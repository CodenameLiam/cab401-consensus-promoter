package qut;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;

public class Test {

    public static void main(String[] args) {
        // Get the current directory for accessing reference genes file
        String path = System.getProperty("user.dir");

        // ----------------------------------------------------------------------------------------------------------//
        // Parallel Program
        // ----------------------------------------------------------------------------------------------------------//
        // Store the start time value
        long parallelStartTime = System.nanoTime();
        // Run the parallel program
        try {
            Parallel.run(path.concat("/referenceGenes.list"), path.concat("/Ecoli"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        // Store the end time
        long parallelEndTime = System.nanoTime();
        // Calculate the total time of execution
        long parallelTotalTime = parallelEndTime - parallelStartTime;
        // Print the time difference
        System.out.println("The parallel program took " + (parallelTotalTime / 1_000_000_000.0) + " seconds to complete.");

        // ----------------------------------------------------------------------------------------------------------//
        // Sequential Program
        // ----------------------------------------------------------------------------------------------------------//
        // Store the start time value
        long sequentialStartTime = System.nanoTime();
        // Run the parallel program
        try {
            Sequential.run(path.concat("/referenceGenes.list"), path.concat("/Ecoli"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        // Store the end time
        long sequentialEndTime = System.nanoTime();
        // Calculate the total time of execution
        long sequentialTotalTime = sequentialEndTime - sequentialStartTime;
        // Print the time difference
        System.out.println("The sequential program took " + (sequentialTotalTime / 1_000_000_000.0) + " seconds to complete.");

        // ----------------------------------------------------------------------------------------------------------//
        // Speedup
        // ----------------------------------------------------------------------------------------------------------//
        double speedUp = sequentialTotalTime / parallelTotalTime;
        System.out.println("The total speedup factor is: " + speedUp);

        // ----------------------------------------------------------------------------------------------------------//
        // Check Results Match
        // ----------------------------------------------------------------------------------------------------------//
        boolean resultsMatch = true;

        HashMap<String, Sigma70Consensus> sequentialResults = Sequential.getConsensus();
        HashMap<String, Sigma70Consensus> parallelResults = Parallel.getConsensus();

        for (String key : sequentialResults.keySet()) {
            String sequentiallKey = sequentialResults.get(key).toString();
            String parallelKey = parallelResults.get(key).toString();
            if (!sequentiallKey.equals(parallelKey)) {
                resultsMatch = false;
                break;
            }
        }

        if (!resultsMatch) {
            System.out.println("The parallel version of the program did not produce the same output as the sequential version");
        } else {
            System.out.println("The parallel and sequential outputs match");
        }
    }
}
