# bioinformatics
to generate the small scripts

### 20220504

Finish the web crawler

Download saccharomyces_cerevisiae genomic fasta files from  http://ftp.ensembl.org/pub/release-106/fasta/saccharomyces_cerevisiae/dna/ using CHECKSUM file 

for i in $(awk '{print $3}' CHECKSUMS); do wget -c  http://ftp.ensembl.org/pub/release-106/fasta/saccharomyces_cerevisiae/dna/$i; done



### 20220505 

Record Everyday Using Unix Script.  Record them in unix.sh

### 20220506
To learn bioawk

To use for loop in shell
