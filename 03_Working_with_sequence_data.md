Nucleotide and protein sequences are stored in two plain-text formats widespread in bioinformatics:**FASTA** and **FASTQ**.

Simple mistakes over minor details like file formats can consume a disproportionate amount of time and energy to discover and fix, so mind these details early on.

#### The FASTA Format
The FASTA format is used to store any sort of sequence data not requiring per-base pair quality scores. This includes reference genome files, protein sequences, coding DNA sequences (CDS), transcript sequences, and so on.

**FASTA** files are composed of sequence entries, each containing two parts: a description and the sequence data. 

#### The FASTQ Format
The FASTQ format extends FASTA by including a numeric quality score to each base in the sequence.
1. The discription line, begining with @
2. Sequence data
3. The line beginning with +, indicates the end of the sequence
4. Quality data, must be the same length as the sequence

#### Nucleotide codes
*Degenerate* nucleotide codes are used to represent two or more bases. Some bioinformatics programsmay handle ambiguous bucleotide differently.

#### Base Quality
1. Be aware of your sequencing technology's error distributions and limitations
2. Consider how this might impact your analyses

#### Inspecting and Trimming Low-Quality Bases
  ## GenomeRanges
  rm(list = ls())
  # Install pacman
  if (!require("pacman", quietly = TRUE))
    install.packages("pacman")
  # Install bioconductor
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  # Install specific Bioconductor Packages
  BiocManager::install(c("qrqc"))
  library(qrqc)
  ######## change the working directory
  #setwd()
  fqfiles <- c(none="untreated1_chr4.fq")
  # Load each file in, using qrqc's readSeqFile
  # We only need qualities, so we turn off some of
  # readSeqFile's other features.
  seq_info <- lapply(fqfiles, function(file) {
    readSeqFile(file, hash=FALSE, kmer=FALSE)
  })
  # Extract the qualities as dataframe, and append
  # a column of which trimmer (or none) was used. This
  # is used in later plots.
  quals <- mapply(function(sfq, name) {
    qs <- getQual(sfq)
    qs$trimmer <- name
    qs
  }, seq_info, names(fqfiles), SIMPLIFY=FALSE)
  # Combine separate dataframes in a list into single dataframe
  d <- do.call(rbind, quals)
  # Visualize qualities
  p1 <- ggplot(d) + geom_line(aes(x=position, y=mean, linetype=trimmer))
  p1 <- p1 + ylab("mean quality (sanger)") + theme_bw()
  print(p1)
  # Use qrqc's qualPlot with list produces panel plots
  # Only shows 10% to 90% quantiles and lowess curve
  p2 <- qualPlot(seq_info, quartile.color=NULL, mean.color=NULL) + theme_bw()
  p2 <- p2 + scale_y_continuous("quality (sanger)")
  print(p2)

#### A FASTA/FASTQ Parsing Example: Counting Nucleotides
This is the beauty of reusing code: wellwritten functions and libraries prevent us from having to rewrite complex parsers.

Reusing software isn’t cheating—it’s how the experts program.

#### Indexed FASTA Files
index the FASTA file
    samtools faidx testfile.fa
access the subsequence for a particular region
    samtools faidx testfile.fa <region>
