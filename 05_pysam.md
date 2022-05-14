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

```
import pysam
samfile = pysam.AlignmentFile("NA12878.dedup.sort.bam","rb")

samfile.check_index()
# return True if index is present

samfile.count(self,contig=None, start=None, stop=None, ...)
# count the number of reads in region
samfile.count('chr1', 10000,10005)
3

samfile.count_coverage(self,contig=None, start=None, stop=None, ...)
# count the number of reads in region, the coverage is computed per-base [ACGT]
# Raises ValueError – if the genomic coordinates are out of range or invalid.
# Returns four array.arrays of the same length in order A C G T
# Return type tuple
samfile.count_coverage('chr1',10000,10005)
(array('L', [0, 2, 2, 0, 0]), array('L', [0, 0, 0, 2, 3]), array('L', [0, 0, 0, 0, 0]), array('L', [1, 0, 0, 0, 0]))

fetch(self, contig=None, start=None, stop=None, ...)
# fetch reads aligned in a region
#Without a contig or region all mapped reads in the file will be fetched. The reads will be returned ordered
#by reference sequence, which will not necessarily be the order within the file. This mode of iteration still
#requires an index. If there is no index, use until_eof=True
samfile.fetch()
<pysam.libcalignmentfile.IteratorRowAllRefs object at 0x101bba978>
samfile.fetch('chr1',10000,10005)
<pysam.libcalignmentfile.IteratorRowRegion object at 0x101bba908>

samfile.find_introns()
# return a dictionary {(start, stop): count} Listing the intronic sites in the reads

get_index_statistics()
# return statistics about mapped/unmapped reads per chromosome as they are stored in the index
samfile.get_index_statistics()
[IndexStats(contig='chrM', mapped=2660, unmapped=0, total=2660), IndexStats(contig='chr1', mapped=179807, unmapped=0, total=179807), IndexStats(contig='chr2', mapped=164021, unmapped=0, total=164021), IndexStats(contig='chr3', mapped=122274, unmapped=0, total=122274), IndexStats(contig='chr4', mapped=101099, unmapped=0, total=101099), IndexStats(contig='chr5', mapped=109653, unmapped=0, total=109653), IndexStats(contig='chr6', mapped=105867, unmapped=0, total=105867), IndexStats(contig='chr7', mapped=111435, unmapped=0, total=111435), IndexStats(contig='chr8', mapped=96166, unmapped=0, total=96166), IndexStats(contig='chr9', mapped=89889, unmapped=0, total=89889), IndexStats(contig='chr10', mapped=103360, unmapped=0, total=103360), IndexStats(contig='chr11', mapped=105766, unmapped=0, total=105766), IndexStats(contig='chr12', mapped=93296, unmapped=0, total=93296), IndexStats(contig='chr13', mapped=54973, unmapped=0, total=54973), IndexStats(contig='chr14', mapped=65766, unmapped=0, total=65766), IndexStats(contig='chr15', mapped=65136, unmapped=0, total=65136), IndexStats(contig='chr16', mapped=84453, unmapped=0, total=84453), IndexStats(contig='chr17', mapped=89884, unmapped=0, total=89884), IndexStats(contig='chr18', mapped=48904, unmapped=0, total=48904), IndexStats(contig='chr19', mapped=81668, unmapped=0, total=81668), IndexStats(contig='chr20', mapped=61529, unmapped=0, total=61529), IndexStats(contig='chr21', mapped=30159, unmapped=0, total=30159), IndexStats(contig='chr22', mapped=49017, unmapped=0, total=49017), IndexStats(contig='chrX', mapped=82231, unmapped=0, total=82231), IndexStats(contig='chrY', mapped=941, unmapped=0, total=941), IndexStats(contig='chrLAM', mapped=0, unmapped=0, total=0)]

```
