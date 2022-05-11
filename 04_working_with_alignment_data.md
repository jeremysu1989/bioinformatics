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

2. 83 # FLAG,the bitwise flag, which contains information abut the alignment. **Bitwise flags**

3. I # RNAME, the reference name, The reference name must be in the SAM/BAM headers as an SQ entry. If the read is unaligned, this entry may be *

4. 2012257 # POS,the position on the reference sequence(using 1-based indexing) of the first mapping base in the query sequence. This may be zero if the reads does not align

5. 40 # MAPQ is the mapping quality, which is a measure of how likely the read is to actually originate from the position it maps to

6. 50M # CIGAR is the *CIGAR* string, which is a specialized format for describing the alignment

7. = # **RNEXT** and **PNEXT** are the reference name and position of a paired-end read's partner.

8. 2011868 # TLEN is the template length for paired-end reads.

-439

9. CAAAAAATTTTGAAAAAAAAAATTGAATAAAAATTCACGGATTTCTGGCT # SEQ stores the original read sequence 

10. 22222222222222222222222222222222222222222222222222 # QUAL stores the origina read base quality

#### Bitwise Flags
Many important pieces of information about an alignment are encoded using *bitwise flags*. Bitwise flags are a very space-efficient and common way to encode attributes, so they're worth understanding. Bitwise flags are much like a series of toggle switches, each of which can be either on or off. Each switch represents whether a particular attribute of an alignment is true or false.
        samtools flags 147

#### CIGAR strings
While bitwise flags store true/false propeties about an alignment, **CIGAR** strings encode information about which bases of an alignment are matches/mismatches, insertions, deletions, soft or hard clipped, and so on.

**Soft clipping** is when only part of the query sequence is aligned to the reference, leaving some portion of the query sequence unaligned. Soft clipping  occurs when an aligner can partially map a read to a location, but the head or tail of the query sequence doesn't match.

**Hard clipping** is similar, but hard-clipped regions are not present in the sequence stored in the SAM field SEQ. A basic CIGAR string contains concatenated pairs of integer *lengths* and character *operations*.  For example, a fully aligned 51 base pair read without insertions or deletions would have a CIGAR string containing a single length/operation pair: **51M**

#### Mapping Qualities
Mapping qualities are one of the most important diagnostics in alignment. All steps downstream of alignment in all bioinformatics projects *critically depend on reliable mapping*.

#### Command-Line Tools for Working with Alignments in the SAM Format
All commands are well documented both online and in the programs themselves.

#### Using **samtools view** to convert between SAM and BAM
Many **samtools** subcommands such as **sort, index, depth,** and **mpileup** all require input files (or streams) to be in BAM format for efficiency.

**samtools view** allows us to convert SAM to BAM with the *-b* option or vice verse
        samtools view -b testfile.sam > testfile.bam
        samtools view testfile.bam > testfile.sam



















