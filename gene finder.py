import argparse

# Function to read the genome from a FASTA file
def read_fasta(file_path):
    with open(file_path, "r") as file:
        genome = ""
        for line in file:
            if not line.startswith(">"):
                genome += line.strip()
    return genome

# Function to find Open Reading Frames (ORFs) in a genome sequence
def find_orfs(genome, start_codon="ATG", stop_codons={"TAA", "TAG", "TGA"}):
    orfs = []
    for frame in range(3):  # For each reading frame
        for i in range(frame, len(genome) - 2, 3):
            codon = genome[i:i+3]
            if codon == start_codon:
                for j in range(i + 3, len(genome) - 2, 3):
                    stop_codon = genome[j:j+3]
                    if stop_codon in stop_codons:
                        orfs.append(genome[i:j+3])
                        break
    return orfs

# Function to generate the reverse complement of a DNA sequence
def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

# Function to find ORFs in both forward and reverse complements (6 frames)
def find_orfs_six_frames(genome):
    orfs = find_orfs(genome)  # Forward frames
    reverse_genome = reverse_complement(genome)
    orfs += find_orfs(reverse_genome)  # Reverse frames
    return orfs

# Function to filter ORFs by length
def filter_by_length(orfs, min_length=100):
    return [orf for orf in orfs if len(orf) // 3 >= min_length]

# Function to check for Shine-Dalgarno sequence (ribosome binding site)
def contains_rbs(sequence, rbs_sequence="AGGAGG", max_distance=20):
    return rbs_sequence in sequence[-max_distance:]

# Function to filter ORFs based on the presence of a Shine-Dalgarno sequence
def filter_by_rbs(orfs, genome, max_distance=20):
    filtered_orfs = []
    for orf in orfs:
        start_index = genome.find(orf)
        if start_index >= max_distance:
            upstream_region = genome[start_index-max_distance:start_index]
            if contains_rbs(upstream_region):
                filtered_orfs.append(orf)
    return filtered_orfs

# Main function to handle command-line arguments and execute the tool
def main():
    parser = argparse.ArgumentParser(description="Gene Finder Tool")
    parser.add_argument("input_file", type=str, help="Path to the FASTA file")
    parser.add_argument("--min_length", type=int, default=100000, help="Minimum length of ORFs (in codons)")
    parser.add_argument("--max_rbs_distance", type=int, default=20, help="Maximum distance for Shine-Dalgarno sequence upstream of start codon")
    args = parser.parse_args()

    genome = read_fasta(args.input_file)
    orfs = find_orfs_six_frames(genome)

    # Filter ORFs by length
    long_orfs = filter_by_length(orfs, args.min_length)

    # Filter ORFs by Shine-Dalgarno sequence (RBS)
    final_orfs = filter_by_rbs(long_orfs, genome, args.max_rbs_distance)

    for orf in final_orfs:
        print(orf)

if __name__ == "__main__":
    main()

#test with a fasta file on powershell, there are 6 files from the NCBIDATASET
