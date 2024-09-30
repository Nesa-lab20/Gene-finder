## README- script/manual WEEK 3-4    
First real git repository — we will write a simple tool that
finds genes in genome sequence using Python, and using a git repository to
track changes as we add features.

# Install and load BioPython on your Terminal

```markdown
ssh baezvg@ilogin.ibex.kaust.edu.sa  
module load bio-python
```

# Create the directory and Repository
The "gene_finder.py" at the end should be one file that works for the whole problem

```markdowm
mkdir Gene-finder-3-4
cd Gene-finder-3-4
git init
touch gene_finder.py README.md
```

# 1 Create the main python file and use an input FASTA file to test
The tool you write needs to take command line parameters
for the input file. The input file should consist of a single FASTA file
containing a single genome.
Your tool should read a FASTA file, and output any region between a
start (‘ATG’) and stop codon (‘TAA’, ‘TAG’, ‘TGA’) in that FASTA file;
you must consider three possible reading frames but may ignore reverse
compliments.

  # Objective: Find Reading Frames from a FASTA file (ecoli.FASTA)
```
nano gene_finder.py

```

# 2 Extend your tool to include the reverse complement and search in all six possible reading frames for genes.


# 3 Use your code to solve the Open Reading Frame problem on Rosalind (Problem 72).
The code was copied to jupyter notebook for the submission on Rosalind
The file with only the code is named and found in this repository as Rosalind72.py
```markdown
Open Rosalind72.py
```

# 4 Simple command for the 14 bacterial genomes- find their ORFs
Apply your code to the genomes you downloaded and find all
Open Reading Frames in the 14 genomes you downloaded.
```
nano extract_orfs.py
python extract_orfs.py
```
Output: The script will create a file named output_proteins.txt containing all the unique open reading frames found in the .fna files in the specified directory.

# 5 Implement a filter by length: discard short ORFs that are unlikely to be functional genes


# 6 Implement gene finder with length, rbs site and rbs type filter
Filter all predicted ORFs based on
whether they contain a Shine-Dalgarno sequence up to 20bp upstream of
the start codon.

