javac -classpath lib/jacobi.jar src/jaligner/matrix/*.java src/jaligner/util/*.java src/jaligner/*.java src/qut/*.java
java -classpath "src;lib/*" qut.Sequential

Open the program in vscode and use Launch.JSON to run

Parallel.Java: 
Parallelised application
Change the parallelisation point (3rd level, 2nd level) in the run method
Change the number of threads using the executor service constructor

Sequential.Java
Sequential version of the application
Run as per normal

Test.Java
Testing area, will run both applications and return the execution time 
Calcualates speedup and prints this to the terminal
If the output does not match, it will return an error message