## Pysam
### Introduction
Pysam is a python module that makes it easy to read and manipulate mapped short read sequence data stored in SAM/BAM files

To use the model to read a file in BAM format, creat a *AlignmentFile* object:
```
>>> import pysam
>>> samfile = pysam.AlignmentFile("NA12878.bam","rb")
```
