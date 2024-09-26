import math

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

# Predicting structure using Chous-Fasman method, (only secondary level)
def chou_fasman(sequence, window_size=6):
    helix_threshold = 1.0
    sheet_threshold = 1.0
    
    structure = []
    
    
    for i in range(len(sequence) - window_size + 1):  # Window size for detecting helices and sheets
        window = sequence[i:i+window_size]
        
        helix_propensity = sum([chou_fasman_data[aa]['helix'] for aa in window]) / window_size
        sheet_propensity = sum([chou_fasman_data[aa]['sheet'] for aa in window]) / window_size
        
        if helix_propensity > helix_threshold:
            structure.append('H')  # Helix
        elif sheet_propensity > sheet_threshold:
            structure.append('S')  # Sheet
        else:
            structure.append('C')  # Coil
    
    return structure

# Smoothing predicted structure because it was effecting number of transitions
def smooth_structure(structure, window_size=3):
    smoothed_structure = []
    for i in range(len(structure)):
        window = structure[max(0, i - window_size // 2) : min(len(structure), i + window_size // 2 + 1)]
        most_common = max(set(window), key=window.count)  # Getting the most common structure in the window
        smoothed_structure.append(most_common)
    return smoothed_structure

# Calculating Shannon entropy
def calculate_shannon_entropy(helix_fraction, sheet_fraction, coil_fraction):
    fractions = [helix_fraction, sheet_fraction, coil_fraction]
    
    # Remove zero probabilities to avoid log(0)
    fractions = [f for f in fractions if f > 0]
    
    entropy = -sum(f * math.log2(f) for f in fractions)
    
    return entropy

# Quantifying Complexity
def quantify_complexity(sequence):
    structure = chou_fasman(sequence)
    
    # Applying smoothing to reduce noise in secondary structure prediction
    smoothed_structure = smooth_structure(structure)
    
    # Variable 1: Secondary Structure Distribution
    helix_count = smoothed_structure.count('H')
    sheet_count = smoothed_structure.count('S')
    coil_count = smoothed_structure.count('C')
    
    total = len(smoothed_structure)
    

    helix_fraction = helix_count / total
    sheet_fraction = sheet_count / total
    coil_fraction = coil_count / total    #important for normalisation
    
    # Shannon entropy calculation based on secondary structure fractions
    entropy = calculate_shannon_entropy(helix_fraction, sheet_fraction, coil_fraction)
    
    # Variable 2: Transitions Between Structures
    transitions = 0
    for i in range(1, len(smoothed_structure)):
        if smoothed_structure[i] != smoothed_structure[i-1]:
            transitions += 1
    
    # Variable 3: Diversity in Secondary Structure Lengths (standard deviation)
    helix_lengths = []
    sheet_lengths = []
    
    current_length = 1
    current_type = smoothed_structure[0]
    
    for i in range(1, len(smoothed_structure)):
        if smoothed_structure[i] == current_type:
            current_length += 1
        else:
            if current_type == 'H':
                helix_lengths.append(current_length)
            elif current_type == 'S':
                sheet_lengths.append(current_length)
            current_type = smoothed_structure[i]
            current_length = 1
    
    if current_type == 'H':
        helix_lengths.append(current_length)
    elif current_type == 'S':
        sheet_lengths.append(current_length)
    
    # Calculating standard deviation for length diversity
    def standard_deviation(lengths):
        if len(lengths) > 1:
            mean = sum(lengths) / len(lengths)
            return math.sqrt(sum((x - mean) ** 2 for x in lengths) / len(lengths))
        return 0
    
    helix_std = standard_deviation(helix_lengths)
    sheet_std = standard_deviation(sheet_lengths)

    # Giving the user options to add weights to each parameter to control results
    entropy_weight = 1.0
    transition_weight = 0.7
    std_weight = 0.1  # Adjusted for standard deviation
    
    complexity_score = (entropy_weight * entropy + 
                        transition_weight * transitions + 
                        std_weight * (helix_std + sheet_std)) / (helix_fraction + sheet_fraction + coil_fraction ) 
                                    #denominator of function acts as a normaliser to counteract
                                     #large values of one variable
    return {
        'helix_fraction': helix_fraction,
        'sheet_fraction': sheet_fraction,
        'coil_fraction': coil_fraction,
        'transitions': transitions,
        'helix_standard_deviation': helix_std,
        'sheet_standard_deviation': sheet_std,
        'shannon_entropy': entropy,
        'complexity_score': complexity_score
    }

# Usage
sequence = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
complexity = quantify_complexity(sequence)

print("Helix Fraction:", complexity['helix_fraction'])
print("Sheet Fraction:", complexity['sheet_fraction'])
print("Coil Fraction:", complexity['coil_fraction'])
print("Number of Transitions:", complexity['transitions'])
print("Helix Standard Deviation:", complexity['helix_standard_deviation'])
print("Sheet Standard Deviation:", complexity['sheet_standard_deviation'])
print("Shannon Entropy:", complexity['shannon_entropy'])
print("Complexity Score:", complexity['complexity_score'])  #returns an arbitrary non standardized value
