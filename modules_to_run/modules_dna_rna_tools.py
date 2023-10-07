from typing import List, Union


def reverse(seqs: str) -> Union[str, List]:
    """
    Reverses the nucleotide sequence
    Arguments:
    - seqs (str): sequence(-s) for processing
    Return:
    - str, for one reversed sequence
    - list, for two or more reversed sequences
    """
    reverse_list = list()
    if len(seqs) <= 1:
        seq_str = seqs[0]
        return seq_str[::-1]
    else:  # if >1 argument is supplied, write the result to list
        for seq in seqs:
            reverse_list.append(seq[::-1])
        return reverse_list


def transcribe(seqs: str) -> Union[str, List]:
    """
    Transcribes the nucleotide sequence
    Arguments:
    - seqs (str): sequence(-s) for processing

    Return:
    - str, for one transcribed sequence
    - list, for two or more transcribed sequences
    """
    transcribe_list = list()
    for seq in seqs:
        seq_str = ''.join(seq)
        for nucleotide in seq_str:
            seq_str = seq_str.replace('T', 'U').replace('t', 'u')
        transcribe_list.append(seq_str)
    if len(transcribe_list) > 1:
        return transcribe_list
    return transcribe_list[0]


def complement(seqs: str) -> Union[str, List]:
    """
    Сomplements the nucleotide sequence
    Arguments:
    - seqs (str): sequence(-s) for processing

    Return:
    - str, for one complemented sequence
    - list, for two or more complemented sequences
    """
    complement_list = list()
    for seq in seqs:
        seq_str = ''.join(seq)
        complement_pairs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A',
                            'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'u': 'a'}
        complement_str = ""
        for nucleotide in seq_str:
            complement_str += complement_pairs[nucleotide]  # go letter by letter and rewrite the corresponding values from the dictionary
        complement_list.append(complement_str)
    if len(complement_list) > 1:
        return complement_list
    return complement_list[0]


def reverse_complement(seqs: str) -> Union[str, List]:
    """
    Reverses and complements the nucleotide sequence
    Arguments:
    - seqs (str): sequence(-s) for processing

    Return:
    - str, for one reversed and complemented sequence
    - list, for two or more reversed and complemented sequences
    """
    if len(seqs) == 1:
        reverse_result = (reverse(seqs))
        reverse_complement_result = ''
        return reverse_complement_result.join(complement(reverse_result))
    else:
        reverse_result_list = reverse(seqs)
        return complement(reverse_result_list)
