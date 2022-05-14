## Pysam
### Introduction
Pysam is a python module that makes it easy to read and manipulate mapped short read sequence data stored in SAM/BAM files

To use the model to read a file in BAM format, creat a ***AlignmentFile*** object:

The bam file should be indexed before import by pysam
```
>>> import pysam
>>> samfile = pysam.AlignmentFile("NA12878.dedup.sort.bam","rb")
```
Once a file is opened you can iterate over all of the reads mapping to a specified region using ***fetch()***. Each iteration returns a ***AlignmentSegment*** object with represents a single read along with its fields and optional tags:
```
>>> for read in samfile.fetch('chr1', 10000, 10005):
...     print(read)
A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT	16	#1	10001	40	147M	*	0	0	TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAA	array('B', [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37])	[('NM', 1), ('MD', '146C0'), ('XM', '...................................................................................................................................................'), ('XR', 'GA'), ('XG', 'GA')]
A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT	0	#1	10002	40	122M	*	0	0	AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCCAACCCTAAA	array('B', [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40])	[('NM', 2), ('MD', '112T8C0'), ('XM', '..HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH..HHH...HHH....'), ('XR', 'CT'), ('XG', 'CT')]
A00265:368:H3M7FDSXY:3:1206:7961:18317:ACTTGAAACGGACTCCTTAC_1:N:0:ACACTAAG+ATCCATAT	0	#1	10005	31	105M	*	0	0	CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCTA	array('B', [40, 40, 37, 40, 37, 40, 40, 40, 40, 40, 37, 40, 40, 37, 37, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 37, 40, 37, 40, 37, 40, 37, 40, 37, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 22, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 37, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 37, 37])	[('NM', 1), ('MD', '103C1'), ('XM', 'HH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHh.'), ('XR', 'CT'), ('XG', 'CT')]
samfile.close()
```

You can also write to a ***AlignmentFile***:
```
import pysam
samfile = pysam.AlignmentFile("ex1.bam", "rb")
reversereads = pysam.AlignmentFile("reversereads.bam", "wb", template=samfile)
  for read in samfile.fetch():
    if read.is_paired:
      reversereads.write(read)
reversereads.close()
samfile.close()
```
