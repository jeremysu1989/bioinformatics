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

#### Finding OVerlapping Ranges

We’ll start with the basic task of finding overlaps between two sets of IRanges objects using the findOverlaps() function

        qry <- IRanges(start=c(1, 26, 19, 11, 21, 7), end=c(16, 30, 19, 15, 24, 8),names=letters[1:6])
        sbj <- IRanges(start=c(1, 19, 10), end=c(5, 29, 16), names=letters[24:26])
        qry
        sbj

![image](https://user-images.githubusercontent.com/104820908/167299704-e389d3c5-de06-4826-8520-332d71ffcdd1.png)

Using the IRanges qry and sbj, we can now find overlaps. Calling findOver laps(qry, sbj) returns an object with class Hits, which stores these overlaps:

        hts <- findOverlaps(qry, sbj)
        hts
![image](https://user-images.githubusercontent.com/104820908/167300016-2acc9848-46d2-4458-a560-f206a1d02ad4.png)

Thinking abstractly, overlaps represent a mapping between query and subject. Depending on how we find overlaps, each query can have many hits in different subjects. A single subject range will always be allowed to have many query hits.

        names(qry)[queryHits(hts)]
        names(sbj)[subjectHits(hts)]
![image](https://user-images.githubusercontent.com/104820908/167300103-299d5dc1-1843-4924-b48e-c2c1dca9f3ec.png)

we could limit our overlap results to only include query ranges that fall entirely within subject ranges with type=within

        hts_within <- findOverlaps(qry, sbj, type="within")
        hts_within

![image](https://user-images.githubusercontent.com/104820908/167300175-48a90001-5aac-4453-9956-a5c937251da9.png)

Another findOverlaps() parameter that we need to consider when computing overlaps is select, which determines how findOverlaps() handles cases where a single query range overlaps more than one subject range. 
Because the options "first", "last", and "arbitrary" all lead findOverlaps() to return only one overlapping subject range per query (or NA if no overlap is found), results are returned in an integer vector where each element corresponds to a query range in qry:

        findOverlaps(qry, sbj, select="first")
        findOverlaps(qry, sbj, select="last")
        findOverlaps(qry, sbj, select="arbitrary")

After running findOverlaps(), we need to work with Hits objects to extract information from the overlapping ranges

1. Hits objects can be coerced to matrix using as.matrix()
2. countQueryHits() returns a vector of how many subject ranges each query IRanges object overlaps. Using the function setNames(),get the resulting vector the same names as original ranges
3. countSubjectHits() and setNames() works as above
4. *create a set of ranges for overlapping regions by calling the ranges() function using the Hits object as the first argument, and the same query and subject ranges we passed to findOverlaps() as the second and third arguments*

        as.matrix(hts)
        countQueryHits(hts)
        setNames(countQueryHits(hts), names(qry))
        countSubjectHits(hts)
        setNames(countSubjectHits(hts), names(sbj))
        ranges(hts, qry, sbj)

The functions subsetByOverlaps() and countOverlaps() simplify some of the most common operations performed on ranges once overlaps are found: keeping only the subset of queries that overlap subjects, and counting overlaps

        countOverlaps(qry, sbj)
        countOverlaps(sbj, qry)
        subsetByOverlaps(qry,sbj)
        subsetByOverlaps(sbj, qry)

#### Finding Nearest Ranges and Calculating Distance
finding ranges that neighbor query ranges

1. nearest() function returns the nearest range, regardless of whether it's upstream or downstream of the query
2. precede() and follow() return the nearest range that the query is upstream or downstream of, respectively

        qry <- IRanges(start=6, end=13, name='query')
        sbj <- IRanges(start=c(2, 4, 18, 19), end=c(4, 5, 21, 24), names=1:4)
        qry
        sbj
        nearest(qry,sbj)
        precede(qry,sbj)
        follow(qry,sbj)

![image](https://user-images.githubusercontent.com/104820908/167326893-7a8e8348-4b8f-46c6-a615-33d71c48a263.png)

these operations are all vectorized

        qry2 <- IRanges(start=c(6, 7), width=3)
        nearest(qry2, sbj)

This family of functions for finding nearest ranges also includes distanceToNearest() and distance(). distanceToNearest() works a lot like findOverlaps() for each query range. distance() returns each pairwise distance between query and subject ranges.

        qry <- IRanges(sample(seq_len(1000), 5), width=50)
        sbj <- IRanges(sample(seq_len(1000), 5), width=50)
        qry
        sbj
        distanceToNearest(qry, sbj)
        distanceToNearest(sbj, qry)
        findOverlaps(qry,sbj)
        distance(qry,sbj)

#### Run Length Encoding and Views
##### Run-length encoding and coverage()

#### Storing Genomic Ranges with GenomicRanges
The GenomicRanges package introduces a new class called GRanges for storing genomic ranges.The GRanges builds off of IRanges. IRanges objects are used to store ranges of genomic regions on a single sequence, and GRanges objects contain the two other pieces of information necessary to specify a genomic location: **sequence name** (e.g., which chromosome) and **strand**.

        library(GenomicRanges)
        gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
                      ranges=IRanges(start=5:8, width=10),
                      strand=c("+", "-", "-", "+"))
        gr

Using the GRanges() constructor, we can also add arbitrary metadata columns by specifying additional named arguments.This illustrates the structure of GRanges objects: **genomic location specified by sequence name, range, and strand (on the left of the dividing bar), and metadata columns (on the right)**. Each row of metadata corresponds to a range on the same row.

        gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
                      ranges=IRanges(start=5:8, width=10),
                      strand=c("+", "-", "-", "+"), gc=round(runif(4), 3))
        gr

We can specif the sequence lengths in the GRanges constructor, or set it after the object has been created using the seqlengths() function.

        seqlens <- c(chr1=152, chr2=432, chr3=903)
        gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
                      ranges=IRanges(start=5:8, width=10),
                      strand=c("+", "-", "-", "+"), gc=round(runif(4), 3),
                      seqlengths=seqlens)
        gr

        seqlengths(gr) <- seqlens
        gr

We access data in GRanges objects much like we access data from IRanges objects: with accessor functions. 

        start(gr) #Ranges accessor function
        end(gr)
        width(gr)
        seqnames(gr) #GRanges specific accessor function
        strand(gr) 
        ranges(gr) # Ranges accessor function
        length(gr)
        names(gr) <- letters[1:length(gr)]
        gr
        table(seqnames(gr))
        gr[seqnames(gr) == "chr1"] # subsetting of GRanges object
        mcols(gr) # access metadata columns, and return a dataframe
        mcols(gr)$gc # access a column of dataframe
        gr$gc

The real power is when we combine subsetting with the data kept in our metadata columns. 

        mcols(gr[seqnames(gr) == "chr1"])$gc
        mean(mcols(gr[seqnames(gr) == "chr1"])$gc)

#### Grouping Data with GRangesList
GRanges objects also have their own version of a list, called GRangesList, which are similar to R’s lists. 

        gr1 <- GRanges(c("chr1", "chr2"), IRanges(start=c(32, 95), width=c(24, 123)))
        gr2 <- GRanges(c("chr8", "chr2"), IRanges(start=c(27, 12), width=c(42, 34)))
        grl <- GRangesList(gr1, gr2)
        grl
        unlist(grl)
        double_grl <- c(grl,grl) # combine many GRangeList objects with c()
        length(double_grl)
        double_grl[1] # returns GRangesList objects
        double_grl[[1]] # returns a list element

GRangesList objects also have some special features. For example, accessor functions for GRanges data (e.g., seqnames(), start(), end(), width(), ranges(), strand(), etc.) also work on GRangesList objects

More often, GRangesLists come about as the result of
using the function split() on GRanges objects. 

        chrs <- c("chr3", "chr1", "chr2", "chr2", "chr3", "chr1")
        gr <- GRanges(chrs, IRanges(sample(1:100, 6, replace=TRUE),
                                    width=sample(3:30, 6, replace=TRUE)))
        seqnames(gr)
        gr_split <- split(gr, seqnames(gr))
        gr_split
        unsplit(gr_split, seqnames(gr))

Grouped data is also the basis of the split-apply-combine pattern 

        lapply(gr_split, function(x) order(width(x)))
        sapply(gr_split, function(x) min(start(x)))
        lapply(gr_split, function(x) min(start(x)))
        sapply(gr_split, length)

for many overlap operation functions (e.g., reduce(), flank(), coverage(), and findOverlaps()), they can work directly with GRangesList objects

#### Working with Annotation Data: GenomicFeatures and rtracklayer

**GenomicFeatures** is a Bioconductor package for creating and working with transcript-based annotation. 

All transcript annotation packages use the same consistent naming scheme—that is, TxDb.<organism>.<annotation-source>.<annotation-version>.

        BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))       
        BiocManager::install(c("TxDb.Mmusculus.UCSC.mm10.ensGene"))  
        library(TxDb.Mmusculus.UCSC.mm10.ensGene)
        txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
        txdb

        #####################################
        TxDb object:
        # Db type: TxDb
        # Supporting package: GenomicFeatures
        # Data source: UCSC
        # Genome: mm10
        # Organism: Mus musculus
        # Taxonomy ID: 10090
        # UCSC Table: ensGene
        # UCSC Track: Ensembl Genes
        # Resource URL: http://genome.ucsc.edu/
        # Type of Gene ID: Ensembl gene ID
        # Full dataset: yes
        # miRBase build ID: NA
        # transcript_nrow: 94647
        # exon_nrow: 348801
        # cds_nrow: 226312
        # Db created by: GenomicFeatures package from Bioconductor
        # Creation time: 2016-09-29 04:15:25 +0000 (Thu, 29 Sep 2016)
        # GenomicFeatures version at creation time: 1.25.17
        # RSQLite version at creation time: 1.0.0
        # DBSCHEMAVERSION: 1.1                







## Bedtools
