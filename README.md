#README- script/manual
#In gene_finder.py, start by setting up a basic main function and include a placeholder for reading command-line arguments:

import argparse

def main():
    parser = argparse.ArgumentParser(description="Gene Finder Tool")
    parser.add_argument("input_file", type=str, help="Path to the FASTA file")
    args = parser.parse_args()

    # Call function to find genes (implementation pending)
    print(f"Reading file: {args.input_file}")

if __name__ == "__main__":
    main()

Commit the initial main functionality

Implement the code to read the FASTA file and find Open Reading Frames (ORFs):

def read_fasta(file_path):
    with open(file_path, "r") as file:
        genome = ""
        for line in file:
            if not line.startswith(">"):
                genome += line.strip()
    return genome

def find_orfs(genome):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    for frame in range(3):
        for i in range(frame, len(genome) - 2, 3):
            codon = genome[i:i+3]
            if codon == start_codon:
                for j in range(i + 3, len(genome) - 2, 3):
                    stop_codon = genome[j:j+3]
                    if stop_codon in stop_codons:
                        orfs.append(genome[i:j+3])
                        break

    return orfs

def main():
    parser = argparse.ArgumentParser(description="Gene Finder Tool")
    parser.add_argument("input_file", type=str, help="Path to the FASTA file")
    args = parser.parse_args()

    genome = read_fasta(args.input_file)
    orfs = find_orfs(genome)

    for orf in orfs:
        print(orf)

if __name__ == "__main__":
    main()

Test your program with a sample FASTA file -6 from the NCBIDATASET bacterial genomes and commit your progress.

Implement reverse complement and ORF search for six frames

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

def find_orfs_six_frames(genome):
    orfs = find_orfs(genome)  # Forward frames
    reverse_genome = reverse_complement(genome)
    orfs += find_orfs(reverse_genome)  # Reverse frames
    return orfs

Extend the Python script to process multiple genomes in one command:

for genome in genomes/*.fasta; do
    python gene_finder.py "$genome"
done

Commit the command for batch processing

Implement ORF Length Filter

def filter_by_length(orfs, min_length=100):
    return [orf for orf in orfs if len(orf) // 3 >= min_length]

Implement Ribosome Binding Site (RBS) Filter

def contains_rbs(sequence, rbs_sequence="AGGAGG", max_distance=20):
    return rbs_sequence in sequence[-max_distance:]

def filter_by_rbs(orfs, genome, max_distance=20):
    filtered_orfs = []
    for orf in orfs:
        start_index = genome.find(orf)
        if start_index >= max_distance:
            upstream_region = genome[start_index-max_distance:start_index]
            if contains_rbs(upstream_region):
                filtered_orfs.append(orf)
    return filtered_orfs

Test your program on multiple FASTA files and genomes to ensure it works as expected

#Command-line arguments: Accepts parameters like input file, minimum ORF length, and RBS search distance.

Command-line Example:

% python gene_finder.py genome.fasta --min_length 150 --max_rbs_distance 15

This will look for ORFs in genome.fasta, filter by ORFs with at least 150 codons, and search for Shine-Dalgarno sequences within 15 bp upstream of the start codon.
