# Gene-finder
# Codon to amino acid table
codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'Stop', 'TAG':'Stop',
    'TGC':'C', 'TGT':'C', 'TGA':'Stop', 'TGG':'W',
}

# Function to compute reverse complement of DNA sequence
def reverse_complement(dna):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join(complement[base] for base in reversed(dna))

# Function to translate a DNA sequence into a protein
def translate_dna(dna):
    protein = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        amino_acid = codon_table.get(codon, '')
        if amino_acid == 'Stop':
            break
        if amino_acid:
            protein.append(amino_acid)
    return ''.join(protein)

# Function to find all distinct ORFs in all six reading frames
def find_orfs(dna):
    orfs = set()
    # Find ORFs in original strand (3 reading frames)
    for frame in range(3):
        for i in range(frame, len(dna)-2, 3):
            if dna[i:i+3] == 'ATG':  # Start codon
                protein = translate_dna(dna[i:])
                if protein:
                    orfs.add(protein)
    
    # Find ORFs in reverse complement strand (3 reading frames)
    reverse_dna = reverse_complement(dna)
    for frame in range(3):
        for i in range(frame, len(reverse_dna)-2, 3):
            if reverse_dna[i:i+3] == 'ATG':  # Start codon
                protein = translate_dna(reverse_dna[i:])
                if protein:
                    orfs.add(protein)
    
    return orfs

# Read input DNA sequence from FASTA format
def read_fasta(data):
    lines = data.strip().split('\n')
    return ''.join(lines[1:])

# Main function to process the input and return the ORFs
def main(fasta_data):
    dna = read_fasta(fasta_data)
    orfs = find_orfs(dna)
    return orfs



