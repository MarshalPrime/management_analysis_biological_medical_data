import math
import csv
import os
from Bio import SeqIO

# Chou-Fasman propensity values
chou_fasman_data = {
    'A': {'helix': 1.42, 'sheet': 0.83, 'coil': 0.66},
    'C': {'helix': 0.70, 'sheet': 1.19, 'coil': 1.19},
    'D': {'helix': 1.01, 'sheet': 0.54, 'coil': 1.46},
    'E': {'helix': 1.51, 'sheet': 0.37, 'coil': 0.74},
    'F': {'helix': 1.13, 'sheet': 1.38, 'coil': 0.60},
    'G': {'helix': 0.57, 'sheet': 0.75, 'coil': 1.56},
    'H': {'helix': 1.00, 'sheet': 0.87, 'coil': 0.95},
    'I': {'helix': 1.08, 'sheet': 1.60, 'coil': 0.47},
    'K': {'helix': 1.16, 'sheet': 0.74, 'coil': 1.01},
    'L': {'helix': 1.21, 'sheet': 1.30, 'coil': 0.59},
    'M': {'helix': 1.45, 'sheet': 1.05, 'coil': 0.60},
    'N': {'helix': 0.67, 'sheet': 0.89, 'coil': 1.56},
    'P': {'helix': 0.57, 'sheet': 0.55, 'coil': 1.52},
    'Q': {'helix': 1.11, 'sheet': 1.10, 'coil': 0.92},
    'R': {'helix': 0.98, 'sheet': 0.93, 'coil': 0.95},
    'S': {'helix': 0.77, 'sheet': 0.75, 'coil': 1.43},
    'T': {'helix': 0.83, 'sheet': 1.19, 'coil': 0.96},
    'V': {'helix': 1.06, 'sheet': 1.70, 'coil': 0.50},
    'W': {'helix': 1.08, 'sheet': 1.37, 'coil': 0.64},
    'Y': {'helix': 0.69, 'sheet': 1.47, 'coil': 1.14}
}

# Chou-Fasman prediction function
def chou_fasman(sequence, window_size=6):
    helix_threshold = 1.0
    sheet_threshold = 1.0
    default_values = {'helix': 0.0, 'sheet': 0.0, 'coil': 0.0}  # Default value for unknown amino acids
    
    structure = []
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        
        helix_propensity = sum([chou_fasman_data.get(aa, default_values)['helix'] for aa in window]) / window_size
        sheet_propensity = sum([chou_fasman_data.get(aa, default_values)['sheet'] for aa in window]) / window_size   #ignore any 'X' units
        
        if helix_propensity > helix_threshold:
            structure.append('H')  # Helix
        elif sheet_propensity > sheet_threshold:
            structure.append('S')  # Sheet
        else:
            structure.append('C')  # Coil
    
    return structure

# Function to smooth the structure predictions with a defined tie-breaking rule
def smooth_structure(structure, window_size=3):
    smoothed_structure = []
    structure_priority = {'H': 1, 'S': 2, 'C': 3}  # Priority order: Helix > Sheet > Coil

    for i in range(len(structure)):
        window = structure[max(0, i - window_size // 2) : min(len(structure), i + window_size // 2 + 1)]
        counts = {aa: window.count(aa) for aa in set(window)}

        # Tie-breaking rule: prioritize H > S > C if counts are the same
        most_common = sorted(counts.items(), key=lambda x: (-x[1], structure_priority[x[0]]))[0][0]
        smoothed_structure.append(most_common)
    
    return smoothed_structure

# Calculate Shannon entropy
def calculate_shannon_entropy(helix_fraction, sheet_fraction, coil_fraction):
    fractions = [helix_fraction, sheet_fraction, coil_fraction]
    fractions = [f for f in fractions if f > 0]
    entropy = -sum(f * math.log2(f) for f in fractions)
    return entropy

# Quantify complexity using Chou-Fasman predictions
def quantify_complexity(sequence):
    structure = chou_fasman(sequence)
    smoothed_structure = smooth_structure(structure)
    
    helix_count = smoothed_structure.count('H')
    sheet_count = smoothed_structure.count('S')
    coil_count = smoothed_structure.count('C')
    total = len(smoothed_structure)

    helix_fraction = helix_count / total
    sheet_fraction = sheet_count / total
    coil_fraction = coil_count / total
    
    entropy = calculate_shannon_entropy(helix_fraction, sheet_fraction, coil_fraction)
    
    transitions = 0
    for i in range(1, len(smoothed_structure)):
        if smoothed_structure[i] != smoothed_structure[i-1]:
            transitions += 1
    
    return {
        'shannon_entropy': entropy,
        'transitions': transitions     #output shannon entropy and transitions
    }

# Function to read sequences from a FASTA file and concatenate them
def read_fasta_sequences(fasta_file):
    concatenated_sequence = ""
    for record in SeqIO.parse(fasta_file, "fasta"):
        concatenated_sequence += str(record.seq)
    return concatenated_sequence

# Function to append the results into an existing CSV file
def append_to_csv(fasta_files, output_csv):
    with open(output_csv, mode='a', newline='') as file:
        writer = csv.writer(file)
        
        for fasta_file in fasta_files:
            concatenated_sequence = read_fasta_sequences(fasta_file)
            print(f"Appended sequence from {fasta_file}: {concatenated_sequence}")
            complexity = quantify_complexity(concatenated_sequence)
            writer.writerow([fasta_file, complexity['shannon_entropy'], complexity['transitions']])

# Main code to list all FASTA files from the folder
def get_fasta_files(folder):
    return [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith('.fasta')]

# Usage
fasta_folder = 'fasta_enzy'  # Directory containing FASTA files
fasta_files = get_fasta_files(fasta_folder)
output_csv = 'output_enzyme.csv'   #output folder

# Check if the file already exists; if not, write the header
try:
    with open(output_csv, 'x', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['FASTA Filename', 'Shannon Entropy', 'Number of Transitions'])
except FileExistsError:
    pass  

append_to_csv(fasta_files, output_csv)

print(f"Results appended to {output_csv}")
