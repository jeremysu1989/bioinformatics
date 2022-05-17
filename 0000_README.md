# bioinformatics
to generate the small scripts

### 20220504

Finish the web crawler

Download saccharomyces_cerevisiae genomic fasta files from  http://ftp.ensembl.org/pub/release-106/fasta/saccharomyces_cerevisiae/dna/ using CHECKSUM file 

for i in $(awk '{print $3}' CHECKSUMS); do wget -c  http://ftp.ensembl.org/pub/release-106/fasta/saccharomyces_cerevisiae/dna/$i; done



### 20220505 

Record Everyday Using Unix Script.  Record them in unix.sh

### 20220506
To learn bioawk ***haven't done*** **Finished at 20220514**

To use for loop in shell 

To check the completeness of all files downloaded from website, use shell script *partial finished, for loop to download dataset from ensemble*  **Finished all** 20220513


### 20220509
To learn **GenomicFeatures** and **rtracklayer**.

### 20220511
To run R script on remote severs

The principle or theory of indexed files in bioinformatics

### 20220513-20220514
Learn pysam

### 20220515
Linux Shell scripting tutorial

### 20220516
Finish bioinformatics data skills
```
import sys
import pysam
from collections import Counter

if len(sys.argv) <2:
  sys.exit("usage: alnstat.py in.bam")
  
fname = sys.argv[1]
bamfile = pysam.AlignmentFile(fname)

stats = Counter()
for read in bamfile:
  stats['total'] += 1
  stats['qcfail'] += int(read.is_qcfail)
  
  # record paired end info
  stats['paired'] += int(read.is_paired)
  stats['read1'] += int(read.is_read1)
  stats['read2'] += int(read.is_read2)
  
  if read.is_unmapped:
    stats['unmapped'] += 1
    continue # other flags dont apply
  
  # recod if mapping quality <=30
  stats['mapping quality <= 30'] += int(read.mapping_quality <= 30)
  
  stats['mapped'] += 1
  stats['proper pair'] += int(read.is_proper_pair)

# specify the output order,since dicts dont have order
output_order = ("total","mapped","unmapped","paired",
                "read1","read2","proper pair","qcfail",
                "mapping quality <= 30")

#format output and print to standard out
for key in output_order:
  format_args = (key, stats[key], 100*stats[key]/float(stats["total"]))
  sys.stdout.write("%s: %d (%0.2f%%)\n" % format_args)

```
