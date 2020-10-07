BOPAL takes two arguments, the newick tree file and the folder containing the sequence files.
As an example, to run BOPAL on the Bacillus data set use the following command:

python main.py tree.dnd Bacillus_Test

To run the pipeline, use the following command (must have OrthoAlign and DupLoss in order to fully run):

python pipeline.py [TESTS FILE] -cne [NUM RUNS] [OUTPUT FOLDER NAME]

TESTS FILE: text file containing set of tests we plan to run (ex. cherry-equal-events-tests.txt).
NUM RUNS: Number of times we want to run each test in the tests file.
OUTPUT FOLDER NAME: Name of output file that we want created that holds all of the results (ex. Test01).

Flags:
-c signals that the tree file given is a cherry
-n signals that the tree file given contains a neighbour (needs to also be a cherry)
-e signals that the events given in the test file have equal probabilities
-d signals that we want to run the DupLoss analysis in our test
-o signals that we want to run the OrthoAlign analysis in our test