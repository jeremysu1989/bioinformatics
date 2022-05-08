## GenomeRanges in R
### A crash course in genomic ranges and coordinate systems
Specify a genomic region or position, we need three necessary pieces of information:
1. Chromosome name
2. Range
3. Strand

ranges are completely linkded to a specific genome version.

Remap genomic range data from an older genomic version's corrdinate system to a newer version's coordinate system
1. CrossMap
2. NCBI Genome Remapping Service
3. LiftOver

Two different range systems
1. 0-based coordinate system, with half-closed,half-open intercals
2. 1-based coordinate system, with closed intervals
![image](https://user-images.githubusercontent.com/104820908/167296071-c6cdf16a-a2d1-43f0-b956-99f6001f59e6.png)



## Bedtools
