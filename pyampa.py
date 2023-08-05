import argparse
import os
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# Torrent index values for each amino acid
torrent_index = {
    'R': 0.106,
    'K': 0.111,
    'C': 0.165,
    'W': 0.172,
    'Y': 0.185,
    'I': 0.198,
    'V': 0.200,
    'H': 0.202,
    'N': 0.240,
    'T': 0.242,
    'F': 0.246,
    'L': 0.246,
    'Q': 0.248,
    'G': 0.265,
    'M': 0.265,
    'S': 0.281,
    'A': 0.307,
    'P': 0.327,
    'E': 0.449,
    'D': 0.479
}

# Add argument parsing
parser = argparse.ArgumentParser(description='Predict antimicrobial peptides.')
parser.add_argument('-i', '--input', required=True, help='Input fasta file')
parser.add_argument('-o', '--output', required=True, help='Output directory')
args = parser.parse_args()

# Read the FASTA file
sequences = list(SeqIO.parse(args.input, 'fasta'))

results = []

# Add a progress bar using tqdm
for seq in tqdm(sequences):
    name = seq.id
    sequence = str(seq.seq)
    
    # Skip sequences that contain non-standard amino acids
    if any(aa not in torrent_index for aa in sequence):
        continue
    
    # Calculate the sliding window index values
    index_values = []
    for i in range(len(sequence) - 6):
        window = sequence[i:i+7]
        window_value = sum(torrent_index[aa] for aa in window) / 7
        if window_value < 0.225:
            index_values.append(1)
        else:
            index_values.append(0)
    
    # Find the longest peptides with only 2 or less negative residues
    i = 0
    while i < len(index_values):
        if sum(index_values[i:i+12]) < 10:  # At least 3 negatives, skip
            i += 1
            continue
        j = i + 12
        while j <= len(index_values) and sum(index_values[i:j]) >= j - i - 2:
            j += 1
        # Here, [i:j-1] is a peptide
        peptide = sequence[i:j-1]
        peptide_average = sum(torrent_index[aa] for aa in peptide) / len(peptide)
        results.append((name, peptide, peptide_average))
        i = j

# Write results to csv
output_directory = args.output
df = pd.DataFrame(results, columns=['Name', 'Peptide', 'Average Torrent Index'])
output_file = os.path.join(output_directory, 'results.csv')
df.to_csv(output_file, index=False)
