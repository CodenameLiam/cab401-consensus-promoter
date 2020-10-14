package qut;

import qut.*;
import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;

public class Parallel {
    private static HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
    // Sigma 70 pattern match will depend on state of the thread
    private static ThreadLocal<Series> sigma70_pattern = ThreadLocal
            .withInitial(() -> Sigma70Definition.getSeriesAll_Unanchored(0.7));
    private static final Matrix BLOSUM_62 = BLOSUM62.Load();
    private static byte[] complement = new byte['z'];

    static {
        complement['C'] = 'G';
        complement['c'] = 'g';
        complement['G'] = 'C';
        complement['g'] = 'c';
        complement['T'] = 'A';
        complement['t'] = 'a';
        complement['A'] = 'T';
        complement['a'] = 't';
    }

    private static List<Gene> ParseReferenceGenes(String referenceFile) throws FileNotFoundException, IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceFile)));
        List<Gene> referenceGenes = new ArrayList<Gene>();
        while (true) {
            String name = reader.readLine();
            if (name == null)
                break;
            String sequence = reader.readLine();
            referenceGenes.add(new Gene(name, 0, 0, sequence));
            consensus.put(name, new Sigma70Consensus());
        }
        consensus.put("all", new Sigma70Consensus());
        reader.close();
        return referenceGenes;
    }

    private static boolean Homologous(PeptideSequence A, PeptideSequence B) {
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f)
                .calculateScore() >= 60;
    }

    private static NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene) {
        int upStreamDistance = 250;
        if (gene.location < upStreamDistance)
            upStreamDistance = gene.location - 1;

        if (gene.strand == 1)
            return new NucleotideSequence(
                    java.util.Arrays.copyOfRange(dna.bytes, gene.location - upStreamDistance - 1, gene.location - 1));
        else {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
            for (int i = 0; i < upStreamDistance; i++)
                result[i] = complement[dna.bytes[reverseStart - i]];

            return new NucleotideSequence(result);
        }
    }

    private static Match PredictPromoter(NucleotideSequence upStreamRegion) {
        return BioPatterns.getBestMatch(sigma70_pattern.get(), upStreamRegion.toString());
    }

    private static void ProcessDir(List<String> list, File dir) {
        if (dir.exists())
            for (File file : dir.listFiles())
                if (file.isDirectory())
                    ProcessDir(list, file);
                else
                    list.add(file.getPath());
    }

    private static List<String> ListGenbankFiles(String dir) {
        List<String> list = new ArrayList<String>();
        ProcessDir(list, new File(dir));
        return list;
    }

    private static GenbankRecord Parse(String file) throws IOException {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }

    public static void run(String referenceFile, String dir) throws FileNotFoundException, IOException {
        // Create a list of callable tasks
        List<Callable<Void>> callableList = new ArrayList<>();

        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);

        for (String filename : ListGenbankFiles(dir)) {
            System.out.println(filename);
            GenbankRecord record = Parse(filename);
            for (Gene referenceGene : referenceGenes) {
                // Parallelising the inner loop
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes) {
                    // Add tasks from the inner for-loop to the callable list
                    callableList.add(new GeneTask(gene, referenceGene, record));
                }

                // // Parallelising the outer loop
                // callableList.add(new referenceTask(record, referenceGene));
            }
        }

        // Runs all of the callable tasks within the callable list
        try {
            executorService.invokeAll(callableList);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet()) {
            System.out.println(entry.getKey() + " " + entry.getValue());
        }
        // Shuts down the executor service
        executorService.shutdown();
    }

    public static void main(String[] args) throws FileNotFoundException, IOException {
        // Get the current directory for accessing reference genes file
        String path = System.getProperty("user.dir");
        // Print action message
        System.out.println("Starting parallel program:");

        // Store the start time value
        long startTime = System.nanoTime();
        // Start the program
        run(path.concat("/referenceGenes.list"), path.concat("/Ecoli"));
        // Store the end time
        long endTime = System.nanoTime();
        // Calculate the total time of execution
        long totalTime = endTime - startTime;
        // Print the time difference
        System.out.println("The parallel program took " + (totalTime / 1_000_000_000.0) + " seconds to complete.");

    }

    // Method to return the consensus results
    public static HashMap<String, Sigma70Consensus> getConsensus() {
        return consensus;
    }

    // Parallelisation Component

    // Create new executor service to handle the assignment of tasks to threads
    // Define thread count in the constructor
    private static ExecutorService executorService = Executors.newFixedThreadPool(24);

    // Add a lock if using this method
    private static final ReentrantLock addLock = new ReentrantLock();

    // Create a callable task to add gene matches to the consensus sequence
    private static class GeneTask implements Callable {

        private Gene gene;
        private Gene referenceGene;
        private GenbankRecord record;

        public GeneTask(Gene gene, Gene referenceGene, GenbankRecord record) {
            this.gene = gene;
            this.referenceGene = referenceGene;
            this.record = record;
        }

        // Critical region using synchronised keyword
        // public synchronized void addMatch(Match prediction) {
        //     consensus.get(referenceGene.name).addMatch(prediction);
        //     consensus.get("all").addMatch(prediction);
        // }

        @Override
        public Object call() {
            if (Homologous(gene.sequence, referenceGene.sequence)) {
                NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                Match prediction = PredictPromoter(upStreamRegion);

                if (prediction != null) {
                    // If using a synchronised method
                    // Synchronise method to ensure flow dependence
                    // addMatch(prediction);

                    // If using a reentrant lock
                    // Lock the thread to ensure flow dependence
                    addLock.lock();
                    consensus.get(referenceGene.name).addMatch(prediction);
                    consensus.get("all").addMatch(prediction);
                    addLock.unlock();
                }
            }
            return null;
        }
    }

    // Create a callable task to parallelise reference
    private static class referenceTask implements Callable {

        private GenbankRecord record;
        private Gene referenceGene;

        public referenceTask(GenbankRecord record, Gene referenceGene) {
            this.record = record;
            this.referenceGene = referenceGene;
        }

        @Override
        public Object call() throws Exception {
            System.out.println(referenceGene.name);

            for (Gene gene : record.genes) {
                // Just parallelising the level 2 loop
                if (Homologous(gene.sequence, referenceGene.sequence)) {
                    NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                    Match prediction = PredictPromoter(upStreamRegion);

                    if (prediction != null) {
                        addLock.lock();
                        consensus.get(referenceGene.name).addMatch(prediction);
                        consensus.get("all").addMatch(prediction);
                        addLock.unlock();
                    }
                }
            }
            return null;
        }
    }
}
