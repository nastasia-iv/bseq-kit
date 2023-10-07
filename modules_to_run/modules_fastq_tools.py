def calculate_average_quality(seq: str) -> int:
    """
    Calculates the average quality of a nucleotide sequence
    Arguments:
    - seq (str): string with values for the quality of each nucleotide

    Return:
    - int, average sequence quality according to the phred33 scale
    """
    converted_sequence_quality = []  # for sequence quality values after conversion
    for base in seq:
        base_quality = (ord(base)) - 33
        converted_sequence_quality.append(base_quality)
    average_quality = (sum(converted_sequence_quality) // len(converted_sequence_quality))
    return average_quality


def calculate_gc_content(seq: str) -> float:
    """
    Calculates the GC composition of a nucleotide sequence
    Arguments:
    - seq (str): sequence to calculate

    Return:
    - int, sequence GC composition
    """
    gc_content = round(((seq.count('G') + seq.count('C')) * 100) / len(seq), 1)
    return gc_content


def calculate_sequence_length(seq: str) -> int:
    """
    Calculates the length of the nucleotide sequence
    Arguments:
    - seq (str): sequence to calculate

    Return:
    - int, sequence length
    """
    len_sequence = len(seq)
    return len_sequence
