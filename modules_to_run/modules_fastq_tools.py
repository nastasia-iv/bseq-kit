import os


def rewrite_fastq_to_dict(input_path: str) -> dict:
    """
     Reads a FASTQ file and returns a dictionary with the reads

     Arguments:
         input_path (str): path to the FASTQ file

     Returns:
         dict: dictionary with the reads in FASTQ format
     """
    seqs = {}  # Инициализируем пустой словарь для хранения данных
    with open(input_path, mode="r") as file:
        for line in file:
            seq_name = line.strip()  # Считываем идентификатор, который является первой строкой
            seq = next(file).strip()  # Считываем последовательность, следующую за идентификатором
            comment = next(file).strip()  # Считываем комментарий, следующий за последовательностью
            seq_quality = next(file).strip()  # Считываем строку с качеством, следующую за комментарием
            seqs[seq_name] = (seq, comment, seq_quality)  # Добавляем кортеж значений в словарь
    return seqs


def calculate_average_quality(seq: str) -> int:
    """
    Calculates the average quality of a nucleotide sequence
    
    Arguments:
        seq (str): string with values for the quality of each nucleotide

    Return:
        int, average sequence quality according to the phred33 scale
    """
    converted_sequence_quality = []  # for sequence quality values after conversion
    for base in seq:
        base_quality = ord(base) - 33
        converted_sequence_quality.append(base_quality)
    average_quality = (sum(converted_sequence_quality) // len(converted_sequence_quality))
    return average_quality


def calculate_gc_content(seq: str) -> float:
    """
    Calculates the GC composition of a nucleotide sequence
    
    Arguments:
        seq (str): sequence to calculate

    Return:
       int, sequence GC composition
    """
    gc_content = round(((seq.upper().count('G') + seq.upper().count('C')) * 100) / len(seq), 1)
    return gc_content


def calculate_sequence_length(seq: str) -> int:
    """
    Calculates the length of the nucleotide sequence
    
    Arguments:
        seq (str): sequence to calculate

    Return:
        int, sequence length
    """
    len_sequence = len(seq)
    return len_sequence


def save_dict_to_fastq(filtered_fastq: dict, input_path: str, output_filename: str) -> str:
    """
     Saves a dictionary with filtered reads to the FASTQ file

     Arguments:
         filtered_fastq (dict): dictionary of data in FASTQ format
         input_path (str): path to the source file
         output_filename (str): file name to save

     Returns:
         str: path where the file was saved
    """
    if output_filename == '':  # Если имя не передано, используем имя входного файла
        output_filename = os.path.basename(input_path)
    if not output_filename.endswith(".fastq"):
        output_filename += ".fastq"
    output_dir = 'fastq_filtrator_results'  # Задаём имя папки для сохранения результата
    os.makedirs(output_dir, exist_ok=True)  # Создаем папку, если она не существует
    output_path = os.path.join(output_dir, output_filename)  # Задаем путь, по к-му сохранится файл
    with open(output_path, mode="w") as filtered_file:  # Запись в файл из словаря
        for seq_id, (sequence, comment, quality) in filtered_fastq.items():  # Запись каждого элемента на отдельную строку
            filtered_file.write(seq_id + '\n')
            filtered_file.write(sequence + '\n')
            filtered_file.write(comment + '\n')
            filtered_file.write(quality + '\n')
    return output_path  # Выводим путь, по к-му сохранился файл
