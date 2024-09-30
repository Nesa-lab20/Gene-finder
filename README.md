## README- script/manual WEEK 3-4    
first real git repository — we will write a simple tool that
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

# Create the main python file and use an input FASTA file to test
The tool you write needs to take command line parameters
for the input file. The input file should consist of a single FASTA file
containing a single genome.
Your tool should read a FASTA file, and output any region between a
start (‘ATG’) and stop codon (‘TAA’, ‘TAG’, ‘TGA’) in that FASTA file;
you must consider three possible reading frames but may ignore reverse
compliments.
```
nano gene_finder.py
```


