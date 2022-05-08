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
![image](https://user-images.githubusercontent.com/104820908/167296147-1179f87a-8cca-4c64-b0fb-ed5bdd9c7951.png)

*You need to mind strand in your work*

### An Interactive Introduction to Range Data with GenomicRanges
#### Installing and Working with Bioconductor Packages
Bioconductor's core packages:
1. GenomicRanges
2. GenomicFeatures
3. Biostrings and BSgenome
4. rtracklayer

        ##GenomeRanges    
        #Install bioconductor
        if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")

        #Install Bioconductor core packages    
        BiocManager::install(version = "3.14")
  
        #Install specific Bioconductor Packages
        BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

#### Storing Generic Ranges with IRanges

        # load the IRanges package
        library(IRanges)

        #create ranges with the IRanges() function
        rng <- IRanges(start = 4, end = 13)
        rng
        IRanges(start=4, width=3)
        IRanges(end=5, width=5)

        # IRanges() constructor (a function that creates a new object) can take vector arguments
        x <- IRanges(start=c(4, 7, 2, 20), end=c(13, 7, 5, 23))
        x

        # each range can be given a name
        names(x) <- letters[1:4]
        x
        class(x)
        str(x)
        start(x)
        end(x)
        width(x)
        end(x) <- end(x) + 4
        names(x)
        range(x)
We can subset IRanges just as we would any other R objects (vectors, dataframes, matrices), using either numeric, logical, or character (name) index.
As with dataframes, indexing using logical vectors created by statements like width(x) > 8 is a powerful way to select the subset of ranges you’re interested in.

        x[2:3]
        start(x) < 5
        x[start(x) < 5]
        x[width(x) > 8]
        x['a']

#### Basic Range Operations: Arithmetic, Transformations, and Set Operations
First, IRanges objects can be grown or shrunk using arithmetic operations like +, -, and *

        x <- IRanges(start=c(40, 80), end=c(67, 114))
        x + 4L
        x - 10L
![image](https://user-images.githubusercontent.com/104820908/167298009-552ef1bf-6ebf-4d94-84c5-64fa060e7b1d.png)

Restrict ranges within a certain bound

        y <- IRanges(start = c(4,6,10,12), width = 13)
        y
        restrict(y,5,10)
![image](https://user-images.githubusercontent.com/104820908/167298111-eb8cab12-3da4-4000-87d5-d1cc01cf9524.png)

flank() is useful in creating ranges upstream and downstream of protein coding genes that could contain promoter sequences

        x
        flank(x, width=7)
        flank(x, width=7, start=FALSE)

![image](https://user-images.githubusercontent.com/104820908/167298535-18c2e9ae-ad78-434b-871b-c25ae7606ee8.png)

reduce() operation takes a set of possibly overlapping ranges and reduces them to a set of nonoverlapping ranges that cover the same positions

        set.seed(0)
        alns <- IRanges(start=sample(seq_len(50), 20), width=5)
        head(alns, 4)
        reduce(alns)
![image](https://user-images.githubusercontent.com/104820908/167298742-caad1a6c-99c5-489a-88c4-0fbd052da219.png)

A similar operation to reduce() is gaps(), which returns the gaps (uncovered portions) between ranges

If you’d like gaps() to include these gaps, specify the start and end positions in gaps (e.g.,gaps(alns, start=1, end=60))

        gaps(alns)
![image](https://user-images.githubusercontent.com/104820908/167298781-a5ad84f5-2d60-46ad-977f-3dc6c33c9f80.png)

Another class of useful range operations are analogous to set operations

        a <- IRanges(start = 4, end = 13)
        b <- IRanges(start = 12, end = 17)
        intersect(a,b)
        setdiff(a,b)
        setdiff(b,a)
        union(a,b)
![image](https://user-images.githubusercontent.com/104820908/167299006-f569fdfc-caaf-42fc-a7b6-f89abd2019ec.png)

IRanges also has a group of set operation functions that act pairwise, taking two equal-length IRanges objects and working range-wise: psetdiff(), pintersect(), punion(), and pgap(). 

















## Bedtools
