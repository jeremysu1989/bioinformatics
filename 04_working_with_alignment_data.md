The most commmon high-throughput data alignmnt format: the sequence alignment/Mapping (**SAM**) format for mapping data (and its binary analog, BAM).

The **SAM** and **BAM** are the standard formats for storing sequencing reads mapped to a reference.

#### Getting to Know Alignment Formats: SAM and BAM

#### The SAM Header
Files in the SAM format consist of **a hearder section** and **an alignment section**.

**Header lines** contain vital metadata about the reference sequences, read and samples information,and (optionally) processing steps and comments.
1. @SQ header entries store information about the reference sequences
2. @RG header entries contain important read group and samples metadata
3. @PG headerentries contain metadata about the programs used to create and process a set of SAM/BAM files

Although *head* works to take a quick peek at the top of a SAM file, kee the following points in mind:
1. *head* won't always provide the entire header
2. it won't work with binary BAM files

The standard way of interacting with SAM/BAM files is through the **SAMtools command-line program (*samtools*)**. 

**samtools view** is the general tool for viewing and converting SAM/BAM files. And usual UNIX tricks can be combinedwith samtools commands.

    samtools view -H testfile.sam
    samtools view -H testfile.bam
    samtools view -H testfile.sam | grep "^@RG"

#### The SAM Alignment Section
Each alignment entry is composed of 11 required fields (and optional fields after this).

    samtools view celegans.sam | tr '\t' '\n' | head -n 11

1. I_2011868_2012306_0:0:0_0:0:0_2489 # QNAME, the query name

2. 83 # FLAG,the bitwise flag, which contains information abut the alignment

3. I # RNAME, the reference name, The reference name must be in the SAM/BAM headers as an SQ entry. If the read is unaligned, this entry may be *

4. 2012257 # POS,the position on the reference sequence(using 1-based indexing) of the first mapping base in the query sequence. This may be zero if the reads does not align

5. 40 # MAPQ is the mapping quality, which is a measure of how likely the read is to actually originate from the position it maps to

6. 50M # CIGAR is the *CIGAR* string, which is a specialized format for describing the alignment

7. = # **RNEXT** and **PNEXT** are the reference name and position of a paired-end read's partner.

8. 2011868 # TLEN is the template length for paired-end reads.

-439

9. CAAAAAATTTTGAAAAAAAAAATTGAATAAAAATTCACGGATTTCTGGCT # SEQ stores the original read sequence 

10. 22222222222222222222222222222222222222222222222222 # QUAL stores the origina read base quality
