#Problem 2 Reverse complements
from Bio import SeqIO
from Bio.Seq import Seq

def find_genes_in_reading_frames(fasta_file):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    
    def find_genes(seq):
        genes = []
        # Check for all reading frames
        for frame in range(3):
            for i in range(frame, len(seq), 3):
                codon = seq[i:i + 3]
                if codon == start_codon:
                    # Start found, search for stop codon
                    for j in range(i, len(seq), 3):
                        stop_codon = seq[j:j + 3]
                        if stop_codon in stop_codons:
                            genes.append(seq[i:j + 3])  # Include the stop codon
                            break
        return genes

    # To hold genes found in all reading frames
    all_genes = []

    # Read sequences from the FASTA file
    sequences = [record.seq for record in SeqIO.parse(fasta_file, "fasta")]

    for seq in sequences:
        # Forward reading frames
        all_genes.extend(find_genes(seq))

        # Reverse complement
        rev_seq = seq.reverse_complement()
        
        # Reverse reading frames
        all_genes.extend(find_genes(rev_seq))

    return all_genes

# Specify your FASTA file path
fasta_file_path = 'extracted_gene.fasta'
genes = find_genes_in_reading_frames(fasta_file_path)

# Print the genes found
for i, gene in enumerate(genes, start=1):
    print(f"Gene {i}:\n{gene}\n")
