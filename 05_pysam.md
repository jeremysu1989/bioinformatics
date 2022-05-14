## Pysam
### 1.1 Introduction
Pysam is a python module that makes it easy to read and manipulate mapped short read sequence data stored in SAM/BAM files

To use the model to read a file in BAM format, creat a ***AlignmentFile*** object:

The bam file should be indexed before import by pysam
```
import pysam
samfile = pysam.AlignmentFile("NA12878.dedup.sort.bam","rb")
```
Once a file is opened you can iterate over all of the reads mapping to a specified region using ***fetch()***. Each iteration returns a ***AlignmentSegment*** object with represents a single read along with its fields and optional tags:
```
for read in samfile.fetch('chr1', 10000, 10005):
    print(read)
A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT	16	#1	10001	40	147M	*	0	0	TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAA	array('B', [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37])	[('NM', 1), ('MD', '146C0'), ('XM', '...................................................................................................................................................'), ('XR', 'GA'), ('XG', 'GA')]
A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT	0	#1	10002	40	122M	*	0	0	AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCCAACCCTAAA	array('B', [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40])	[('NM', 2), ('MD', '112T8C0'), ('XM', '..HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH..HHH...HHH....'), ('XR', 'CT'), ('XG', 'CT')]
A00265:368:H3M7FDSXY:3:1206:7961:18317:ACTTGAAACGGACTCCTTAC_1:N:0:ACACTAAG+ATCCATAT	0	#1	10005	31	105M	*	0	0	CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCTA	array('B', [40, 40, 37, 40, 37, 40, 40, 40, 40, 40, 37, 40, 40, 37, 37, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 37, 40, 37, 40, 37, 40, 37, 40, 37, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 22, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 37, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 37, 37])	[('NM', 1), ('MD', '103C1'), ('XM', 'HH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHh.'), ('XR', 'CT'), ('XG', 'CT')]
samfile.close()
```

You can also write to a ***AlignmentFile***:
```
import pysam
samfile = pysam.AlignmentFile("NA12878.dedup.sort.bam","rb")
reversereads = pysam.AlignmentFile("reversereads.bam", "wb", template=samfile)
for read in samfile.fetch():
    if read.is_paired:
        reversereads.write(read)
reversereads.close()
samfile.close()
```

An alternative way of accessing the data in a SAM file is by iterating over each base of a specified region using the ***pileup()*** method. Each iteration returns a ***PileupColumn*** which represents all the reads in the SAM file that mao to a single base in the reference sequence. This list of reads are represented as ***PileupRead*** objects in the ***PileupColumn.pileups*** property:
```
import pysam
samfile = pysam.AlignmentFile("NA12878.dedup.sort.bam","rb")
for pileupcolumn in samfile.pileup('chr1', 10000, 10005):
    print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:                                     
            print('\tbase in read %s = %s' % (pileupread.alignment.query_name,pileupread.alignment.query_sequence[pileupread.query_position]))
samfile.close()

# The above code outputs:

coverage at base 10000 = 1
	base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = T

coverage at base 10001 = 2
	base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = A
	base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = A

coverage at base 10002 = 2
	base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = A
	base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = A

coverage at base 10003 = 2
	base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C
	base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C

coverage at base 10004 = 3
	base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C
	base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C
	base in read A00265:368:H3M7FDSXY:3:1206:7961:18317:ACTTGAAACGGACTCCTTAC_1:N:0:ACACTAAG+ATCCATAT = C

coverage at base 10005 = 3
	base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C
	base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C
	base in read A00265:368:H3M7FDSXY:3:1206:7961:18317:ACTTGAAACGGACTCCTTAC_1:N:0:ACACTAAG+ATCCATAT = C
...
```

Commonds available in *samtools* are available as simple function calls.
```
pysam.sort("-o", "output.bam", "ex1.bam")
# corresponds to the commond line:
samtools sort -o output.bam ex1.bam
```

## 1.2 API (Application Programming Interface)
### 1.2.1 SAM/BAM/CRAM files
Objects of type *AlignmentFile* allow working with BAM/SAM formatted files.

**Class** pysam.**AlignmentFile**
	AlignmentFile(filepath_or_object, mode=None. template=None, reference_name=None, ...)
