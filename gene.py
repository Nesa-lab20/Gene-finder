module load bio-python

nano gene-finder.py 


from Bio import SeqIO

# Replace 'ecoli.fna' with the correct path to your file if necessary
for record in SeqIO.parse("ecoli.fna", "fasta"):
    print("ID:", record.id)
    print("Sequence:", record.seq)
    print("Description:", record.description)

from Bio import SeqIO
import glob

# Use glob to find all .fna files in the specified directory
for file_path in glob.glob("ncbi_dataset/data/*.fna"):
    for record in SeqIO.parse(file_path, "fasta"):
        print("File:", file_path)
        print("ID:", record.id)
        print("Sequence:", record.seq)
        print("Description:", record.description)
        print()  # Print a newline for better readability
----
