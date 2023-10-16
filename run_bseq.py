from typing import List, Union
from modules_to_run.modules_dna_rna_tools import is_dna_rna, reverse, transcribe, complement, reverse_complement
from modules_to_run.modules_amino_acid_tools import is_peptide, operations_with_aminoacid
from modules_to_run.modules_fastq_tools import rewrite_fastq_to_dict, calculate_average_quality, calculate_gc_content, calculate_sequence_length, save_dict_to_fastq


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
        is_dna_rna(seq)
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


def run_amino_acid_tools(*seqs: str, operation: str) -> list:
    """
    Run AminoAcid Tools (to calculate the molecular weight of a sequence(-s) or find out the percentage of amino acids)

    Arguments:
    - *seqs (str): one or more string sequences to be analyzed
    - operation (str): action to be done with sequence(-s)

    Return:
    - list with the result of operation
    """
    if operation not in operations_with_aminoacid:
        raise ValueError('Incorrect operation')
    output = []
    for seq in seqs:
        is_peptide(seq)
        output.append(operations_with_aminoacid[operation](seq))
    return output


def run_fastq_tools(input_path: str, output_filename: str = '', gc_bounds: Union[tuple, int] = (0, 100), length_bounds: Union[tuple, int, float] = (0, 2**32), quality_threshold: int = 0) -> str:
    """
    Run FastqTools (to checks each fastq sequence in the dict 'seqs' against the specified conditions)

    Arguments:
    - input_path (str): path to file with fastq sequences to check
    - output_filename (str): optional, name of the file to save the result
    - gc_bounds (Union[tuple, int, float]): GC composition interval (percents) to filter, default is (0, 100); if one value is passed, it's considered the upper limit
    - length_bounds (Union[tuple, int]): sequence length interval to filter, default is (0, 2**32); if one value is passed, it's considered the upper limit
    - quality_threshold (int): threshold value of average read quality (phred33) to filter, default is 0

    Return:
    - str: path for filtered file
    """
    # Обрабатываем границы ГЦ-состава
    if isinstance(gc_bounds, (float, int)):
        gc_min = 0
        gc_max = gc_bounds
    else:
        gc_min = gc_bounds[0]
        gc_max = gc_bounds[1]
    # Обрабатываем границы по длине
    if isinstance(length_bounds, int):
        length_min = 0
        length_max = length_bounds
    else:
        length_min = length_bounds[0]
        length_max = length_bounds[1]
    
    seqs = rewrite_fastq_to_dict(input_path)  # переводим fastq в словарь

    filtered_reads = {}

    for seq_name, values in seqs.items():
        seq, comment, seq_quality = values

        gc_content = calculate_gc_content(seq)
        length = calculate_sequence_length(seq)
        avg_quality = calculate_average_quality(seq_quality)

        if (
                gc_min <= gc_content <= gc_max
                and length_min <= length <= length_max
                and avg_quality >= quality_threshold
        ):
            filtered_reads[seq_name] = values

    return save_dict_to_fastq(filtered_reads, input_path, output_filename)
