#This code was generated with the assistance of ChatGPT, powered by OpenAI's GPT-4 architecture. The model is designed to understand and generate human-like text based on the input it receives. Below is the full prompt that was used to guide the development of the provided Python code for extracting distinct protein strings from a DNA sequence in FASTA format:

#Prompt: "Create a Python script that reads a DNA sequence from a given FASTA format input. The script should extract all distinct candidate protein strings that can be translated from open reading frames (ORFs) of the DNA sequence. Ensure the script includes a codon table for translation and handles both the original and reverse complement strands. The output should be a sorted list of distinct protein strings, displayed in a specified order. Provide the full code along with a sample input and expected output."

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
    return DNA_CODON_TABLE.get(codon) #Returns Nonen if codon is not found

def reverse_complement(dna):
    lookup = str.maketrans('ATGC', 'TACG') #translation table for complement
    return dna.translate(lookup)[::-1] #Reverse and complement in one step

def possible_protein_strings(s):
    results = set()
    length = len(s)

    for i in range(length - 2):
        codon = s[i:i + 3]
        if translate_codon(codon) == 'M':  # Start translation at 'M'
            protein_string = ''
            for j in range(i, length - 2, 3):
                protein = translate_codon(s[j:j + 3])
                if protein == 'Stop':
                    if protein_string:  # Ensure we have a valid protein before adding
                        results.add(protein_string)
                    break
                if protein:
                    protein_string += protein

    return results

def read_fasta(data):
    lines = data.strip().split('\n')
    return ''.join(lines[1:])

if __name__ == "__main__":
    # Input in FASTA format
    fasta_data = """>Rosalind_8967
CTCCGTTTTAAACACTTAGACGAAACACCAACTTGATCTCTGTTAGTCTGAATGTGTTTC
GTCCGTACCGCGCCCTCCTGTACAACGCAAGTGTATGTTTAAGTTGAACGAGTCGAAACA
ACTCTGTTGCTCCTAAGGGTTCGCCGAAACGGCCCGACTGTGAGGGGGACGAGGAGTAGG
TTCTAACTCTGCAATCTTGTCCATATTCTGCGGAGATTGCATGAAAACCAAGACAAAAGG
TATTATCTTCTAGGACCCGAAGGCCGGCCGTCCTGGCGGCTACGGTGAACTGAAGGCTCA
ATACCAGCAGAGGACGCTCTATGCTATAGTAGCGCGAGCGATCGTCTACGTCTGTAGTGG
GTTGGTACGTCTGCTGGGTTGATCTGGTGTATCACTATAGTGAAGCCTGCAACTGTCGGG
AGGTGCGCAGACGCGAATGAATTGCACTGCGAAGGCACTGTAGCTACAGTGCCTTCGCAG
TGCAATTCATTCAAGCACTTACCCTCGACTGTAAACTTATGTTCGCCGGTACTGGAAACG
CCAGACGATTATCCAGATTATGGAGTCTCCGATTACCACAGAGCTGTACCACAAAAGTGA
ACAGGGTCGTTGGTGATGGCGACGATGGGTATGCGGCGTAACCTGCGTGAGCGCTCCCAA
TACCATTGAGGTGAGGAGGCTTCCGAGCTGGAAGTGTGTAGTGGCGACAACTGAAAAACC
TAGTCATCCATGCAGTGGGTTGAATGGAAACAGGATGGCTCGGAATCGGTCCTTTAGACA
ATTGGAGTCCGTGCGATAGATTGAATCTGATTATGGAGTTATGTACGATTGTCGCCGTGT
CACCTAAACTGCACTTTATCGAATAACTGCGGTGCTGCGCCCTTTCCTTGGAACCCCAAC
ATCACCTACTCCTCTCTATGTAGTTG"""
    
    # Read the DNA sequence
    dna_sequence = read_fasta(fasta_data)
    
    # Generate reverse complement
    reverse_dna_sequence = reverse_complement(dna_sequence)
    
    # Combine proteins from both strands
    proteins = possible_protein_strings(dna_sequence) | possible_protein_strings(reverse_dna_sequence)
    
    # Sort the results to ensure the correct order(this not required)
    sorted_proteins = sorted(proteins, key=lambda x: (len(x), x))

    # Print the results
    for protein in sorted_proteins:
        print(protein)
