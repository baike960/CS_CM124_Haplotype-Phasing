In this project, the goal is to develop and implement an algorithm that takes as an input the genotypes of individuals coming from admixed populations and infers the haplotypes of the individuals.
Our algorithm performs Expectation Maximization (EM) on a sequence of
several SNPs for each individual at a time, returns a phase with maximum probability, and moves on to the next SNP sequence of the same window size until all SNPs are covered. We use parsimony approach on an overlap region of each two consecutive sequences to arrange and connect the chunk phases.


1. To run the code and generate solution files for test data, put all four files under the same directory as the test data. 
2. Open chunk_ligation_v2.py, go to the second last line and change the string value inside construct_genotype as the name of the test data. Change the name of the solution file generated on the last line if needed. Save it and exit. 
2. On Unix terminal, use the command line "python2 chunk_ligation_v2.py" to run the code. 
