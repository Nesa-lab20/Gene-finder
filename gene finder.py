
#main code
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
#test with a fasta file on powershell, there are 6 files from the NCBIDATASET
