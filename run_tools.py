from typing import List, Union

from modules_to_run.modules_dna_rna_tools import reverse, transcribe, complement, reverse_complement
from modules_to_run.modules_amino_acid_tools import operation_dict
from modules_to_run.modules_fastq_tools import calculate_average_quality, calculate_gc_content, calculate_sequence_length

nucleotides = {'A', 'T', 'G', 'C', 'U', 'a', 't', 'g', 'c', 'u'}


def is_dna(seq: str) -> bool:
    """
    Checks whether the string is DNA/RNA or not

        Arguments:
    - seq (str): sequence(-s) for processing

        Return:
    - bool, if string is DNA/RNA; throws an error if it not
    """
    t_count = 0
    u_count = 0
    for nucleotide in seq:
        if nucleotide not in nucleotides:
            raise ValueError('Incorrect nucleotide sequence')
        if nucleotide.upper() == 'T':
            t_count += 1
        if nucleotide.upper() == 'U':
            u_count += 1
    if t_count != 0 and u_count != 0:
        raise ValueError('Incorrect nucleotide sequence')
    return True


def run_dna_rna_tools(*seqs: str, operation: str = '') -> Union[str, List]:
    """
    Run DnaRnaTools (to performs standard operations with nucleotide sequences)

    Arguments:
    - seqs (str): one or more sequences for processing
    - operation (str): name of the operation to be performed (reverse, transcribe, complement, reverse_complement)

    Return:
    - str, for one sequence to be processing
    - list, for two or more sequences to be processing
    """
    for seq in seqs:
        is_dna(seq)
    # start the required operation
    if operation == 'reverse':
        return reverse(seqs)
    if operation == 'transcribe':
        return transcribe(seqs)
    if operation == 'complement':
        return complement(seqs)
    if operation == 'reverse_complement':
        return reverse_complement(seqs)
    raise ValueError('Incorrect operation')


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
    if set(seq).issubset(one_letter_aminoacids) is True:
        return True
    raise ValueError('Incorrect amino acid sequence')


def run_amino_acid_tools(*seqs: str, operation: str = '') -> list:
    """
    Run AminoAcid Tools (to calculate the molecular weight of a sequence(-s) or find out the percentage of amino acids)

    Arguments:
    - *seqs (str): one or more string sequences to be analyzed
    - operation (str): action to be done with sequence(-s)

    Return:
    - list with the result of operation
    """
    if operation not in operation_dict:
        raise ValueError('Incorrect operation')
    output = []
    for seq in seqs:
        is_peptide(seq)
        output.append(operation_dict[operation](seq))
    return output


def run_fastq_tools(seqs: dict, gc_bounds: Union[tuple, int] = (0, 100), length_bounds: Union[tuple, int, float] = (0, 2**32), quality_threshold: int = 0) -> dict:
    """
    Run FastqTools (to checks each fastq sequence in the dict 'seqs' against the specified conditions)

    Arguments:
    - seqs (dict): fastq sequences to check
    - gc_bounds (Union[tuple, int, float]): GC composition interval (percents) to filter, default is (0, 100); if one value is passed, it's considered the upper limit
    - length_bounds (Union[tuple, int]): sequence length interval to filter, default is (0, 2**32); if one value is passed, it's considered the upper limit
    - quality_threshold (int): threshold value of average read quality (phred33) to filter, default is 0

    Return:
    - dict, including fastq sequences that have passed all the conditions
    """
    # process gc_bounds conditions
    if isinstance(gc_bounds, (float, int)):
        gc_min = 0
        gc_max = gc_bounds
    else:
        gc_min = gc_bounds[0]
        gc_max = gc_bounds[1]
    # process length_bounds conditions
    if isinstance(length_bounds, int):
        length_min = 0
        length_max = length_bounds
    else:
        length_min = length_bounds[0]
        length_max = length_bounds[1]

    gc_filtered_seqs: dict = {}
    length_filtered_seqs: dict = {}
    all_filtered_seqs: dict = {}
    for seq_name, seq_and_quality in seqs.items():
        if gc_min <= calculate_gc_content(seqs.get(seq_name)[0]) <= gc_max:
            gc_filtered_seqs[seq_name] = seq_and_quality
    for seq_name, seq_and_quality in gc_filtered_seqs.items():
        if length_min <= calculate_sequence_length(seqs.get(seq_name)[0]) <= length_max:
            length_filtered_seqs[seq_name] = seq_and_quality
    for seq_name, seq_and_quality in length_filtered_seqs.items():
        if quality_threshold <= calculate_average_quality(seqs.get(seq_name)[-1]):
            all_filtered_seqs[seq_name] = seq_and_quality
    return all_filtered_seqs
