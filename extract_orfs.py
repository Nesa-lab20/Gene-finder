import os

DNA_CODON_TABLE = {
    'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',
    'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
    'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',
    'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
    'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
    'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
    'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
    'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
    'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',
    'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
    'TAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
    'TAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
    'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',
    'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
    'TGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
    'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
}

def translate_codon(codon):
    return DNA_CODON_TABLE.get(codon)

def reverse_complement(dna):
    lookup = str.maketrans('ATGC', 'TACG')
    return dna.translate(lookup)[::-1]

def possible_protein_strings(s):
    results = set()
    length = len(s)

    for i in range(length - 2):
        codon = s[i:i + 3]
        if translate_codon(codon) == 'M':
            protein_string = ''
            for j in range(i, length - 2, 3):
                protein = translate_codon(s[j:j + 3])
                if protein == 'Stop':
                    if protein_string:
                        results.add(protein_string)
                    break
                if protein:
                    protein_string += protein

    return results

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.read().strip().split('\n')
        return ''.join(lines[1:])  # Skip the header

def process_files(directory):
    all_proteins = set()
    
    for filename in os.listdir(directory):
        if filename.endswith('.fna'):
            file_path = os.path.join(directory, filename)
            dna_sequence = read_fasta(file_path)
            reverse_dna_sequence = reverse_complement(dna_sequence)
            
            proteins = possible_protein_strings(dna_sequence) | possible_protein_strings(reverse_dna_sequence)
            all_proteins.update(proteins)

    return all_proteins

if __name__ == "__main__":
    directory = 'ncbi_dataset/data'
    proteins = process_files(directory)
    
    with open('output_proteins.txt', 'w') as output_file:
        for protein in sorted(proteins):
            output_file.write(protein + '\n')
