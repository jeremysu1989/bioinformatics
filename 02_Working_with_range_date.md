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











## Bedtools
