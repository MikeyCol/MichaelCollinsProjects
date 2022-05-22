These files cotain command line programs that have various methods 
of analysing tRNA sequences when input in the FASTA format.

bos-tRNA.fa 
this is the test input for findUnique.py
contains 22 tRNA sequences in FASTA format

tass2.fa
this is the test input for findORFs.py
contains one long DNA sequence

findUnique.py
This finds all unique continuous substrings of each of the 22 tRNAs and prints them to stdout

sequenceAnalysis.py
Has various tools for analysing DNA/RNA sequences including the orfFinder function. 
orfFinder finds the protein created by the RNA sequence on each open reading frame (orf)

findORFs.py
Comandlline shell for running orfFinder

compare.py
checks two text files for inconsistencies and prints the first pair of inconsistent lines
and the line number at the inconsitency. 



