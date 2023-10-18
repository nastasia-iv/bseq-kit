from typing import Dict

one_letter_aminoacids = {
    'A', 'R', 'N', 'D', 'C', 'H', 'G', 'Q', 'E', 'I',
    'L', 'K', 'M', 'P', 'S', 'Y', 'T', 'W', 'F', 'V'
}


def is_peptide(seq: str) -> bool:
    """
    Check whether the incoming sequence is an aminoacid
    Arguments:
        - seq (str): sequence for processing

    Return:
        - bool: result of the check
    """
    return set(seq).issubset(one_letter_aminoacids)


def calculate_amino_acid_percentage(seq: str) -> Dict[str, float]:
    """
    Calculates the percentage of amino acids in the entered amino acid sequence

    Arguments:
    	- seq (str): amino acid sequences to be analyzed

    Return:
    	- str: a string with the percentage of each amino acid
    """
    amino_acid_counts: Dict[str, int] = {}  # dict to store count of each amino acid
    for amino_acid in seq:
        if amino_acid in amino_acid_counts:
            amino_acid_counts[amino_acid] += 1
        else:
            amino_acid_counts[amino_acid] = 1
    total_amino_acids = len(seq)
    amino_acid_percentages = {}  # dict to store each amino acid and its %
    for amino_acid, count in amino_acid_counts.items():
        percentage = round(((count / total_amino_acids) * 100), 2)
        amino_acid_percentages[amino_acid] = percentage
    return amino_acid_percentages


def calculate_molecular_weight(seq: str) -> float:
    """
    Calculates the molecular weight of entered amino acid sequence

    Arguments:
    	- seq (str): amino acid sequences to be analyzed
    Return:
    	- float: molecular weight value for amino acid sequence
    """
    amino_acid_weights = {
        'G': 57.051, 'A': 71.078, 'S': 87.077, 'P': 97.115, 'V': 99.131,
        'T': 101.104, 'C': 103.143, 'I': 113.158, 'L': 113.158, 'N': 114.103,
        'D': 115.087, 'Q': 128.129, 'K': 128.172, 'E': 129.114, 'M': 131.196,
        'H': 137.139, 'F': 147.174, 'R': 156.186, 'Y': 163.173, 'W': 186.210
    }
    weight = 18.02  # for the H and OH at the termini
    for amino_acid in seq:
        weight += amino_acid_weights[amino_acid]
    return round(weight, 2)


def calculate_hydrophobicity_eisenberg(seq: str) -> float:
    """
    Calculate estimation of hydrophilicity/hydrophobicity of amino acid sequence
    
    Arguments:
    	- seq (str): amino acid sequences to be analyzed
    	
    Return:
    	- float: rough estimation of hydrophilicity/hydrophobicity of amino acid sequence
    """
    hydrophobicity_values = {
        'A': 0.5, 'R': 0.65, 'N': 1.0, 'D': 1.3, 'C': -0.15,
        'Q': 1.0, 'E': 1.5, 'G': 0.75, 'H': 0.7, 'I': -1.3,
        'L': -1.3, 'K': 0.75, 'M': -1.1, 'F': -1.9, 'P': 0.55,
        'S': 0.6, 'T': 0.3, 'W': -0.5, 'Y': -1.65, 'V': -0.9
    }
    hydrophobicity_sum = sum(hydrophobicity_values.get(amino_acid, 0) for amino_acid in seq)
    return round(hydrophobicity_sum, 2)


def calculate_pI(seq: str) -> float:
    """
    Calculate the isoelectric point of aminoacid sequence
    
    Arguments:
    	- seq (str): amino acid sequences to be analyzed
    	
    Return:
    	- float: isoelectric point of amino acid sequence
    """
    pK_values = {
        'A': (2.34, 9.60),
        'R': (2.17, 9.04, 12.48),
        'N': (2.02, 8.80),
        'D': (2.09, 9.82, 3.86),
        'C': (1.71, 8.33, 10.30),
        'Q': (2.17, 9.13),
        'E': (2.19, 9.76, 4.25),
        'G': (2.34, 9.60),
        'H': (1.82, 9.17, 6.00),
        'I': (2.32, 9.76),
        'L': (2.36, 9.60),
        'K': (2.18, 8.95, 10.5),
        'M': (2.28, 9.21),
        'F': (2.58, 9.24),
        'P': (2.00, 10.60),
        'S': (2.21, 9.15),
        'T': (2.63, 10.43),
        'W': (1.22, 9.39),
        'Y': (2.20, 9.11, 10.10),
        'V': (2.29, 9.72)
    }
    # Initialization of variables for leftmost and rightmost elements
    N_end_pK = 0
    C_end_pK = 0
    # Find the marginal elements and their corresponding pKs
    for amino_acid in seq:
        if amino_acid in pK_values:
            pK_list = pK_values[amino_acid]
            if len(pK_list) >= 2:
                if N_end_pK == 0:
                    N_end_pK = pK_list[1]  # Второй pK
                C_end_pK = pK_list[0]  # Первый pK
    # Calculate pI
    total_pK = N_end_pK + C_end_pK
    count = 2  # We take into account the found pKs - there are at least 2
    for amino_acid in seq:
        if amino_acid in pK_values:
            pK_list = pK_values[amino_acid]
            if len(pK_list) >= 3:
                total_pK += pK_list[2]  # Третий pK
                count += 1
    pI = total_pK / count
    return round(pI, 2)


operations_with_aminoacid = {
    'calculate_molecular_weight': calculate_molecular_weight,
    'calculate_amino_acid_percentage': calculate_amino_acid_percentage,
    'calculate_hydrophobicity_eisenberg': calculate_hydrophobicity_eisenberg,
    'calculate_pI': calculate_pI
}
