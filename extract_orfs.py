
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
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith('.fna'):
                file_path = os.path.join(dirpath, filename)
                dna_sequence = read_fasta(file_path)
                print(f"Processing file: {file_path}")  # Debug: show full path
                print(f"DNA Sequence: {dna_sequence}")  # Debug: show sequence
                
                reverse_dna_sequence = reverse_complement(dna_sequence)
                proteins = possible_protein_strings(dna_sequence) | possible_protein_strings(reverse_dna_sequence)
                print(f"Found proteins: {proteins}")  # Debug: show found proteins
                
                # Create an output file for the current FASTA file
                output_file_name = os.path.splitext(filename)[0] + '_output_proteins.txt'
                output_file_path = os.path.join(dirpath, output_file_name)
                
                with open(output_file_path, 'w') as output_file:
                    for protein in sorted(proteins):
                        output_file.write(protein + '\n')

                if not proteins:
                    print("No proteins found in this file.")

if __name__ == "__main__":
    directory = 'ncbi_dataset/data'
    process_files(directory)
