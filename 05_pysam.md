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

samfile.get_reference_length('chr1')
# return reference length corresponding to numerical tid
249250621

samfile.get_reference_name(24)
# return reference name corresponding to numerical tid
'chrY'
samfile.get_reference_name(23)
'chrX'
samfile.get_reference_name(22)
'chr22'
samfile.get_reference_name(0)
'chrM'

samfile.get_tid('chr1')
1

samfile.getrname(2)
'chr2'

samfile.gettid('chr4')
4

samfile.has_index()
# 类似于这样的函数可以用在python中进行条件判断而决定后续脚本是否执行
True

samfile.head(1)
# return an iterator over the first n alignments
<pysam.libcalignmentfile.IteratorRowHead object at 0x101bbaac8>
samfile.head(2)
<pysam.libcalignmentfile.IteratorRowHead object at 0x101bbaa58>

samfile.is_valid_reference_name('chr1')
True

samfile.is_valid_tid(3)
True

samfile.lengths
(16571, 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 48502)

samfile.references
('chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrLAM')

samfile.mapped
# int with total number of mapped alignments according to the statistics recorded in the index
2099954

samfile.nocoordinate
# int with total number of reads without coordinates according to the statistics recorded in the index
0

samfile.nreferences
# int with the number of reference sequences in the file
26

samfile.pileup('chr1',10000,10005)
# perform a pileup within a region
<pysam.libcalignmentfile.IteratorColumnRegion object at 0x101abf480>
```

**Class** pysam.**AlignmentHeader**

header information for a AlignmentFile object

**Class** pysam.**AlignmentSegment**

class representing an aligned segment

**Class** pysam.**AlignmentPileupColumn**

A pileup of reads at a particular reference sequence position (column). A pileup column contains all the reads that map to a certain target base.

### 1.2.2 Tabix files

*TabixFile* opens tabular files that have been indexed with **tabix**

**Class** pysam.**TabixFile**

Random access to bgzf formatted files that have been indexed by tabix

### 1.2.3 FASTA files

**Class** pysam.**FastaFile**

Random access to fasta formatted files that have been indexed by faidx

### 1.2.4 FASTQ files

**Class** pysam.**FastxFile**

Stream access to fasta or fastq formatted files

```
# two different way to load a file
fastq="NA12878_R0.fq"
Reads = pysam.FastxFile(fastq)
for id in Reads:
    print(id.name)
# another way to load
with pysam.FastxFile("NA12878_R0.fq") as fh:
    for entry in fh:
    	print(entry.name)
```

extract the umi and index info from fastq
```
with pysam.FastxFile("NA12878_R1.fq") as R2, open("umi_index.txt", mode='w') as fout:
    for entry in R2:
        fout.write(str(entry.comment) + '\n')
	# or
	# fout.write(str(entry.comment) + "UMI_index" + str(len(entry.sequence)) + '\n')
```

### 1.2.5 VCF/BCF files

### 1.2.6 HTSFile

## 1.3 Working with BAM/CRAM/SAM-formatted files
### 1.3.1 Opening a file

### 1.3.2 Fetching reads mapped to a region
Reads are obtained through a call to the pysam.AlignmentFile.fetch() method which returns an iterator.
Each call to the iterator will returns a pysam.AlignedSegment object
```
samfile = pysam.AlignmentFile("NA12878.dedup.sort.bam")
iter = samfile.fetch("chr1", 10000,10005)
for x in iter:
    print(x)

A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT	16	#1	10001	40	147M	*	0	0	TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAA	array('B', [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37])	[('NM', 1), ('MD', '146C0'), ('XM', '...................................................................................................................................................'), ('XR', 'GA'), ('XG', 'GA')]
A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT	0	#1	10002	40	122M	*	0	0	AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCCAACCCTAAA	array('B', [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40])	[('NM', 2), ('MD', '112T8C0'), ('XM', '..HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH..HHH...HHH....'), ('XR', 'CT'), ('XG', 'CT')]
A00265:368:H3M7FDSXY:3:1206:7961:18317:ACTTGAAACGGACTCCTTAC_1:N:0:ACACTAAG+ATCCATAT	0	#1	10005	31	105M	*	0	0	CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCTA	array('B', [40, 40, 37, 40, 37, 40, 40, 40, 40, 40, 37, 40, 40, 37, 37, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 37, 40, 37, 40, 37, 40, 37, 40, 37, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 22, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 37, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 40, 40, 37, 40, 40, 40, 40, 40, 37, 37])	[('NM', 1), ('MD', '103C1'), ('XM', 'HH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHH...HHh.'), ('XR', 'CT'), ('XG', 'CT')]
```
### 1.3.3 Using the pileip-engine
In contrast to **fetching**, the **pileup** engine returns for each base in the **reference** sequence the reads that map to that particular position. In the typical view of reads stacking vertically on top of the reference sequence similar to a multiple alignment, **fetching** iterates over the rows of this implied multiple alignment while a **pileup** iterates over the **columns**.

Calling **pileup()** will return an iterator over each **column** (reference base) of a specified region. 
```
import pysam
samfile = pysam.AlignmentFile("NA12878.dedup.sort.bam","rb")
outfile = open("base_position.txt", mode='w')
for read in samfile.pileup('chr1', 10000, 10005):
     outfile.write("\ncoverge at base %s = %s" % (read.pos, read.n))
     for base in read.pileups:
          outfile.write("\tbase in read %s = %s" % (base.alignment.query_name,base.alignment.query_sequence[base.query_position]))
outfile.close()
samfile.close()

less base_position.txt 
coverge at base 10000 = 1       base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = T
coverge at base 10001 = 2       base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = A  base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = A
coverge at base 10002 = 2       base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = A  base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = A
coverge at base 10003 = 2       base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C  base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C
coverge at base 10004 = 3       base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C  base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C  base in read A00265:368:H3M7FDSXY:3:1206:7961:18317:ACTTGAAACGGACTCCTTAC_1:N:0:ACACTAAG+ATCCATAT = C
coverge at base 10005 = 3       base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C  base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = C  base in read A00265:368:H3M7FDSXY:3:1206:7961:18317:ACTTGAAACGGACTCCTTAC_1:N:0:ACACTAAG+ATCCATAT = C
coverge at base 10006 = 3       base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = T  base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = T  base in read A00265:368:H3M7FDSXY:3:1206:7961:18317:ACTTGAAACGGACTCCTTAC_1:N:0:ACACTAAG+ATCCATAT = T
coverge at base 10007 = 3       base in read A00265:368:H3M7FDSXY:3:2569:7093:19914:TACGGGTAAGCGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = A  base in read A00265:368:H3M7FDSXY:3:2475:13991:6464:AAGCTGTCTCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT = A  base in read A00265:368:H3M7FDSXY:3:1206:7961:18317:ACTTGAAACGGACTCCTTAC_1:N:0:ACACTAAG+ATCCATAT = A
...
```

### 1.3.4 Creating BAM/CRAM/SAM files from scratch
```
import pysam
header = { 'HD': {'VN': '1.0'},
           'SQ': [{'LN': 1575, 'SN': 'chr1'},
                  {'LN': 1584, 'SN': 'chr2'}] }
with pysam.AlignmentFile("created_test.bam", "wb", header=header) as outf:
     a = pysam.AlignedSegment()
     a.query_name = "read_28833_29006_6945"
     a.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
     a.flag = 99
     a.reference_id = 0
     a.reference_start = 32
     a.mapping_quality = 20
     a.cigar = ((0,10), (2,1), (0,25))
     a.next_reference_id = 0
     a.next_reference_start=199
     a.template_length=167
     a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
     a.tags = (("NM", 1),
               ("RG", "L1"))
     outf.write(a)

samtools view created_test.bam 
read_28833_29006_6945	99	chr1	33	20	10M1D25M	=	200	167	AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG	<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<	NM:i:1	RG:Z:L1
```

### 1.3.5 Using streams

### 1.3.6 Using samtools commands within python

Commands available in samtools are available as simple function calls. Command line options are provided as arguments
```
pysam.sort("-o", "output.bam", "ex1.bam")
# corresponds to the command line
samtools sort -o output.bam ex1.bam
```
