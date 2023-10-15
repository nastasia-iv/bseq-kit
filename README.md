# bseq-kit

## Overview 

***bseq-kit*** is a mini-program that allows the user to work with nucleotide and amino acid sequences, as well as fastq format data. This tool contains four semantic and functional blocks:

#### :round_pushpin: run_dna_rna_tools
Block for working with nucleotide sequences. It can converts a DNA/RNA sequence into its reverse, complement, reverse-complement form, or translate to its trancribe form.

#### :round_pushpin: run_amino_acid_tools
Block for working with amino acid sequences. Calculates the molecular weight of a sequence or its amino acid composition.

#### :round_pushpin: run_fastq_tools
Block for working with fastq format. Filters fastq-sequences based on user-specified conditions.

#### :round_pushpin: bio_files_processor
Block for searching sequences in a gbk file and for converting multi-line fasta files into single-line ones.

## Usage

To use *run_** modules simply import `run_bseq` into your `*.py` script as shown below:
```python
import run_bseq
```


To use `bio_files_processor` simply import it into your `*.py` script as shown below:
```python
import bio_files_processor
```

### :inbox_tray: Input

#### run_dna_rna_tools

The program has two required input parameters:
* Nucleotide sequence(s), variable argument. Sequences to be analyzed are **not** case-specific (`str` type);
* `operation = 'option'` Name of the operation to be executed, keyword argument (`str` type):

```python
run_bseq.run_dna_rna_tools('AGTtC', 'GTT', operation = 'reverse')  # correct
```

:exclamation: You must use one of predefined operation names described in the "Options" section below.

#### run_amino_acid_tools

The program has two required input parameters:
* Aminoacid sequence(s), variable argument. The amino acid sequence must be written in single-letter form and in uppercase (`str` type);
* `operation = 'option'` Name of the operation to be executed, keyword argument (`str` type):

```python
run_bseq.run_amino_acid_tools('ARDF', operation = 'calculate_molecular_weight')  # correct
```

:exclamation: You must use one of predefined operation names described in the "Options" section below.


#### run_fastq_tools

The program has one required input parameter:
* `input_path` Fastq-file with reads that need to be filtered  (`str` type). 

The remaining arguments indicate filtering conditions. For filtering parameters you can leave the default values or change it.  Operation to filtering described in the "Options" section below.

```python
run_fastq_tools(input_path = '/work_dir/example.fastq', gc_bounds = (20, 50), length_bounds = (90), quality_threshold = 30)  # correct
```

#### bio_files_processor

##### convert_multiline_fasta_to_oneline
* `input_fasta` 
Path to fasta-file for convert (`str` type). 

* `output_fasta`
Optional argument. A new fasta file in which each sequence fits on one line. If the `output_fasta` parameter is not passed, then the name of the input file is taken instead. If you do not add a file extension to `output_fasta`, the function will automatically add the extension `.fasta`.

##### select_genes_from_gbk_to_fasta
* `input_gbk`
Path to the input GBK file.

* `genes` 
Genes of interest, next to which neighbors are searched.
:exclamation:  At least 2 genes must be passed on.

* `n_before`, `n_after` 
Number of genes before and after (>0). Default values are 1.

* `output_fasta` 
Name of the output file. Optional argument. 
If the `output_fasta` parameter is not passed, then the name of the input file is taken instead. If you do not add a file extension to `output_fasta`, the function will automatically add the extension `.fasta`.


### :outbox_tray: Output

#### run_dna_rna_tools

`str` If one sequence is given as input, `list` if two or more sequences are submitted:

```python
run_bseq.run_dna_rna_tools('AGTtC', 'GTT', operation = 'reverse')
['CtTGA', 'TTG']

run_bseq.run_dna_rna_tools('GGTca', operation = 'transcribe')
'GGUca'
```

#### run_amino_acid_tools

`list` For both one and several arguments. 
*Please note*, when using the `calculate_amino_acid_percentage` operation, a `list[dict]` is always returned as output:

```python
run_bseq.run_amino_acid_tools('ARDF', operation = 'calculate_amino_acid_percentage')
[{'A': 25.0, 'R': 25.0, 'D': 25.0, 'F': 25.0}]

run_bseq.run_amino_acid_tools('ARDF', operation = 'calculate_molecular_weight')  
[507.54]
```

#### run_fastq_tools

New fastq-file which consisting only of those sequences that passed all the specified conditions. The filtered data will be saved to the `fastq_filtrator_results` folder (if it does not exist, it will be created automatically). 

 #### bio_files_processor
 
* #### convert_multiline_fasta_to_oneline
Fasta file in which each sequence fits into one line.

 * #### select_genes_from_gbk_to_fasta
Fasta file containing protein sequences of gene neighbors for specified genes of interest.



## Options

Each program block contains its own implementable operations.


#### run_dna_rna_tools  

 -   `transcribe`  
 
Convert nucleotide (DNA/RNA) sequence into its transcribe form.
*Please note*, entering an RNA sequence to run this operation will produce the original sequence.

-   `reverse`  

Convert nucleotide (DNA/RNA) sequence into reverse sequence.

-   `complement`  

Convert nucleotide (DNA/RNA) sequence into complement sequence.

-   `reverse_complement`  

Convert nucleotide (DNA/RNA) sequence into complement sequence and reverse it.

#### run_amino_acid_tools  

-   `calculate_molecular_weight`

Calculate the molecular weight of the input amino acid sequences based on the mass of each amino acid residue. Reference values for the masses of amino acid residues are taken from [The University of Washington's Proteomics Resource (UWPR)](https://proteomicsresource.washington.edu/protocols06/masses.php) and rounded to three decimal places. The calculations took into account the presence of *H* and *OH* groups at the termini of the sequences. The output is molecular weight of sequence in Daltons (the result is rounded to two decimal places).

-   `calculate_amino_acid_percentage`  

Calculate the percentage of each amino acid in the sequence. The input is a string with amino acid sequence. The output is the percentage of each amino acid in the sequence (the result is rounded to two decimal places).

- `calculate_pI`

Calculate the isoelectric point of amino acid sequences. The function operation is based on the formula for determining the isoelectric point:

$$pI = \dfrac{(pK_1 + pK_2 + ... + pK_n)}{n},$$

where $pK$ is dissociation constant of free $NH_2$ and $COOH$ radicals in amino acids.
The output is the value of the isoelectric point of aminoacids sequence (the result is rounded to two decimal places).

- `calculate_hydrophobicity_eisenberg`

Calculate estimation of hydrophilicity/hydrophobicity of amino acid sequences. The function operation is based on the Einzenberg hydrophilicity/hydrophobicity scale of amino acids. The output is the rough estimate of the hydrophobicity of aminoacids sequence (the result is rounded to two decimal places). A value less than 0 indicates hydrophobicity; a value greater than 0 indicates hydrophilicity.

#### run_fastq_tools  

-   `output_filename`  

Optional argument. Name for the new file with filtered reads.
If the `output_filename` parameter is not passed, then the name of the input file is taken instead. If you do not add a file extension to `output_filename`, the function will automatically add the extension `.fastq`.

-   `gc_bounds`  

Composition GC interval (in percent) for filtering. Default is `tuple` (0, 100). When passing a single number as an argument, it is considered an upper bound. For example, `gc_bounds = (20, 80)` save only reads with GC content from 20 to 80%, `gc_bounds = 50`  save reads with GC content less than or equal to 50%.
:exclamation: When specifying only the upper bound, use `int` or `float`

-   `length_bounds` 

Length interval for filtering.  Default is `tuple` (0, 2**32). When passing a single number as an argument, it is considered an upper bound. For example, `length_bounds = (150, 500)` save only reads with length from 150 to 500 bp, `length_bounds = 400`  save reads with length less than or equal to 400 bp.
:exclamation: When specifying only the upper bound, use `int` or `float`

-   `quality_threshold` 

The threshold value of average read quality for filtering. Default is 0 (`int`) in Phred33 scale. Reads with average quality across all nucleotides below the threshold will be discarded.


## Troubleshooting

For `run_dna_rna_tools`, the program automatically checks whether the input is nucleotide sequence (contains only letters *A*, *T*, *G*, *C*, *U* in any case and not *T* and *U* together). 
For `run_amino_acid_tools` the correctness of the amino acid sequence is checked (whether it contains only letters corresponding to one-letter amino acid abbreviations).
The name of the operation is also checked. If an error is detected, the program will be interrupted and a corresponding message will be displayed on the screen:

```python
run_bseq.run_dna_rna_tools('AGTUU', operation = 'reverse')  # incorrect input sequence
ValueError: Incorrect sequence

run_bseq.run_amino_acid_tools('ARDF', operation = 'percentage')  # incorrect operation
ValueError: Incorrect operation
```
#


