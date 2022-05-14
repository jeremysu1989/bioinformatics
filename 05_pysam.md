## Pysam
### Introduction
Pysam is a python module that makes it easy to read and manipulate mapped short read sequence data stored in SAM/BAM files

To use the model to read a file in BAM format, creat a ***AlignmentFile*** object:

The bam file should be indexed before import by pysam
```
>>> import pysam
>>> samfile = pysam.AlignmentFile("NA12878.bam","rb")
```
Once a file is opened you can iterate over all of the reads mapping to a specified region using ***fetch()***. Each iteration returns a ***AlignmentSegment*** object with represents a single read along with its fields and optional tags:
```
```
