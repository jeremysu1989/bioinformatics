The most commmon high-throughput data alignmnt format: the sequence alignment/Mapping (**SAM**) format for mapping data (and its binary analog, BAM).

The **SAM** and **BAM** are the standard formats for storing sequencing reads mapped to a reference.

#### Getting to Know Alignment Formats: SAM and BAM

##### The SAM Header
Files in the SAM format consist of **a hearder section** and **an alignment section**.

**Header lines** contain vital metadata about the reference sequences, read and samples information,and (optionally) processing steps and comments.
1. @SQ header entries store information about the reference sequences
2. @RG header entries contain important read group and samples metadata
3. @PG headerentries contain metadata about the programs used to create and process a set of SAM/BAM files

Although *head* works to take a quick peek at the top of a SAM file, kee the following points in mind:
1. *head* won't always provide the entire header
2. it won't work with binary BAM files

The standard way of interacting with SAM/BAM files is through the **SAMtools command-line program (*samtools*)**. 
