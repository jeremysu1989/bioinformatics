Nucleotide and protein sequences are stored in two plain-text formats widespread in bioinformatics:**FASTA** and **FASTQ**.

Simple mistakes over minor details like file formats can consume a disproportionate amount of time and energy to discover and fix, so mind these details early on.

#### The FASTA Format
The FASTA format is used to store any sort of sequence data not requiring per-base pair quality scores. This includes reference genome files, protein sequences, coding DNA sequences (CDS), transcript sequences, and so on.

**FASTA** files are composed of sequence entries, each containing two parts: a description and the sequence data. 

#### The FASTQ Format
The FASTQ format extends FASTA by including a numeric quality score to each base in the sequence.
1. The discription line, begining with @
2. Sequence data
3. The line beginning with +, indicates the end of the sequence
4. Quality data, must be the same length as the sequence
