from Bio import SeqIO

def find_reading_frames(fasta_file):
    # Read sequences from the FASTA file
    sequences = [record.seq for record in SeqIO.parse(fasta_file, "fasta")]

    reading_frames = []

    for seq in sequences:
        # Forward frames
        for frame in range(3):
            reading_frames.append(seq[frame:].translate())

        # Reverse complement the sequence
        rev_seq = seq.reverse_complement()
        
        # Reverse frames
        for frame in range(3):
            reading_frames.append(rev_seq[frame:].translate())

    return reading_frames

# Specify your FASTA file path
fasta_file_path = 'extracted_gene.fasta'
frames = find_reading_frames(fasta_file_path)

# Print the reading frames
for i, frame in enumerate(frames, start=1):
    print(f"Reading Frame {i}:\n{frame}\n")

