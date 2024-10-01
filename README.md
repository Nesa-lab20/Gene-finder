## README- script/manual WEEK 3-4    
```markdown
#LLM used: CHAT GPT-40 mini
Below is the full prompt that was used to guide the development of the provided Python code for extracting distinct protein strings from a DNA sequence in FASTA format:
#Prompt: 
#"Create a Python script that reads a DNA sequence from a given FASTA format input. The script should extract all distinct candidate protein strings that can be translated from open reading frames (ORFs) of the DNA sequence. Ensure the script includes a codon table for translation and handles both the original and reverse complement strands"
```
## START 
First real git repository — we will write a simple tool that
finds genes in genome sequence using Python, and using a git repository to
track changes as we add features.

>**NOTE= Most of the .txt files for the output were too big to clone on to Github**

# Install and load BioPython on your Terminal

```markdown
ssh baezvg@ilogin.ibex.kaust.edu.sa  
module load bio-python
```

# Create the directory and Repository
The "gene_finder.py" at the end should be one file that works for the whole problem- ANSWER: the complete file SDS.py is the end result

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

  # Objective: Find Reading Frames from a FASTA file (extracted_gene.FASTA)
```
nano gene_finder.py
python gene_finder.py

```
The file with the answer for applying the code gene_finder.py to the extracted_gene.FASTA is ORFS.txt

# 2 Extend your tool to include the reverse complement and search in all six possible reading frames for genes.
We can just modify the code gene_finder.py
Copy the code and make a new script

```
nano gene_finder2.py
python gene_finder2.py

```
The file with the answer for applying the code gene_finder.py to the ecoli.FASTA is ORF2.txt

# 3 Use your code to solve the Open Reading Frame problem on Rosalind (Problem 72).
The LLM was used to modify the code "gene_finder2.py" to accomodate the input given by Rosalind.
The code was copied to jupyter notebook for the submission on Rosalind
The file with only the code is named and found in this repository as Rosalind72.py
```markdown
Open Rosalind72.py
```

# 4 Simple command for the 14 bacterial genomes- find their ORFs
Apply your renewed code from Rosalind to the genomes you downloaded from the ncbidataset and find all
Open Reading Frames in the 14 genomes you downloaded.
You need to know where are the .fna files in the directory
```
nano extract_orfs.py
python extract_orfs.py
```
Output: The script will create a file named output_proteins.txt containing all the unique open reading frames found in the .fna files in the specified directory.

# 5 Implement a filter by length: discard short ORFs that are unlikely to be functional genes
Modify your code to minimize the ORFs found, the less likely to be functional genes

```
nano DiscardShortORFS.py
python DiscardShortORFS.py
```
Output: The script will create a file named filtered_proteins_min_length.txt containing all the unique open reading frames, without the discarded short ones.

# 6 Implement gene finder with length, rbs site and rbs type filter
Filter all predicted ORFs based on
whether they contain a Shine-Dalgarno sequence up to 20bp upstream of
the start codon.
```
nano SDS.py
python SDS.py
```
Output: The script will create a file named filtered_proteins_with_SD.txt containing all the unique open reading frames found in the .fna files in the specified directory.

