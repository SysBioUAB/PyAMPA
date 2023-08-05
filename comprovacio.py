from Bio.SeqUtils.ProtParam import ProteinAnalysis
import math

def calc_features(sequence):
    # Define the set of nonpolar residues
    nonpolar_residues = set('ACGILMFPWYV')

    # Calculate the percentage of nonpolar residues
    NP = sum(1 for aa in sequence if aa in nonpolar_residues) / len(sequence)
    print(NP)

    # Calculate the presence of Trp and Tyr
    W = 1 if 'W' in sequence else 0
    Y = 1 if sequence.count('Y') >= 2 else 0

    # Calculate the isoelectric point
    protein_analysis = ProteinAnalysis(sequence)
    IP_val = protein_analysis.isoelectric_point()
    IP = 1 if IP_val > 10 else 0

    return NP, W, Y, IP

def calc_half_life(sequence):
    # Calculate the peptide features
    NP, W, Y, IP = calc_features(sequence)

    # Calculate the natural log of the half-life
    ln_t_half = 2.226 + (0.053 * NP) - (1.515 * W) + (1.290 * Y) - (1.052 * IP)

    # Convert to the half-life
    t_half = math.exp(ln_t_half)

    return t_half


sequence = "ATQLGHKLGRKKK"
pI = calc_half_life(sequence)
print(f"Half-life of the sequence '{sequence}': {pI}")