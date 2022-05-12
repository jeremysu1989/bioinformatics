## Analysing bisulfite methylation sequencing data

### Bisulfite-Seq theory and Quality Control

#### Distribution of CG
- CpG dinuleotides are not evenly distributed
- Most occurwithin high density regions called CpG Islands
- Most of the genome contains methylated CpGs
- CpGs in CpG islands are largely unmethylated

#### Regulation by DNA methylation
Faults in correct DNA methylation may result in
- early development failure
- epigenetic syndromes
- cancer

#### Bisulphite Library Preparation
- Whole Genome Bisulphite Sequencing
- Post Bisulphite Adapter Tagging
- Reduced Representation (RRBS)

#### Bisulfite conversion of a genomic locus
- 2 different PCR products
- 4 possible different sequence strands from one genomics locus
- each of these 4 sequence strands can theoretically exist in any possibles conversion state

#### Commond sequencing protocols
- Directional libraries **OT OB**
- PBAT libraries **CTOT CTOB**
- Non-directional libraries **OT OB CTOT CTOB**

#### BS-seq ANalysis Workflow
![image](https://user-images.githubusercontent.com/104820908/168039076-dd98c6a0-d9cd-4945-9955-303ee542926b.png)

#### Part I : Intitial QC - What does QC tell you about your library?
- #of sequences
- Basecall qualities
- Base composition
- Potential contaminants
- Expected duplication rate

**Imortant to trim because failre to do so might result in**
``` diff
- Low mapping efficiency
+ Mis-alignments
! Errros in methylation cals since adapters are methylated
# Basecall errors tend toward 50%(C:mC)
```
#### PartII : Sequence alignment - Bismark primary alignment output (BAM file)
- chromosome
- position
- sequence
- quality
- methylation call

#### Sequence duplicaiton
![image](https://user-images.githubusercontent.com/104820908/168040635-372105e8-0f71-43a4-a79e-034fc6f0ddbb.png)

#### Part III : Mapped QC - Methylation bias
**M-Bias Plot** good opportunity to look at conversion efficiency

#### Bismark workflow
**Pre Alignment**

FastQC

Trim Galore

**Alignment**

Bismark

**Post Alignment**

Deduplication

Methylation extractor

bismark2report

#### Useful links

![image](https://user-images.githubusercontent.com/104820908/168083820-653504ac-5f31-435f-9f89-abd9ce41d432.png)
---

## Exercises: QC and Mapping of BS-Seq data
**Software**
- FastQC
- Trim Galore
- Cutadapt
- Bismark
- Bowtie2

#### Exerciese outline
![image](https://user-images.githubusercontent.com/104820908/168051685-d9f0fa4e-8158-4f5f-96ca-9b27ec694bf7.png)

#### Deduplication
```
deduplicate_bismark NA12878-20200324-L01.bam
```
#### Methylation extraction
```
bismark_methylation_extractor --bedGraph --gzip NA12878-20200324-L01.deduplicated.sam
```

This commond
- generates strand and context specific cytosine output files
- counts overlapping parts of read pairs only once
- generates an M-bias report
- produces a bedGraph and coverage file
- generates an overall count report for (splitting report)

**Output files**
```diff
CpG_OT_NA12878-20200324-L01.deduplicated.txt.gz
CpG_OB_NA12878-20200324-L01.deduplicated.txt.gz
CpG_CTOT_NA12878-20200324-L01.deduplicated.txt.gz
CpG_CTOB_NA12878-20200324-L01.deduplicated.txt.gz
CHH_OT_NA12878-20200324-L01.deduplicated.txt.gz
CHH_OB_NA12878-20200324-L01.deduplicated.txt.gz
CHH_CTOT_NA12878-20200324-L01.deduplicated.txt.gz
CHH_CTOB_NA12878-20200324-L01.deduplicated.txt.gz
CHG_OT_NA12878-20200324-L01.deduplicated.txt.gz
CHG_OB_NA12878-20200324-L01.deduplicated.txt.gz
CHG_CTOT_NA12878-20200324-L01.deduplicated.txt.gz
CHG_CTOB_NA12878-20200324-L01.deduplicated.txt.gz
NA12878-20200324-L01.deduplicated_splitting_report.txt
NA12878-20200324-L01.deduplicated.M-bias.txt
NA12878-20200324-L01.deduplicated.bismark.cov.gz
NA12878-20200324-L01.deduplicated.bedGraph.gz
```

#### Bismark HTML report
```
bismark2report
```

**Output file**
```
NA12878-20200324-L01_SE_report.html
```

---

## Visualising and Exploring BS-Seq Data

```
samtools view NA12878-20200324-L01.bam | wc -l
23106639

samtools view NA12878.dedup.bam | wc -l
2099954

wc -l C*
49862706 total

```

**Reads summary for Bismark Methylation Extraction**
| file_name        | Total line No | Uniq line No |
|-----------------:|--------------:|-------------:|
| CHG_CTOB_NA12878| 6914286| 984351| 
| CHG_CTOT_NA12878| 7003633|   998226| 
| CHG_OB_NA12878| 9684|     8048| 
| CHG_OT_NA12878| 10372|     8680| 
| CHH_CTOB_NA12878| 15055662|   991894| 
| CHH_CTOT_NA12878| 15252260|  1005774| 
| CHH_OB_NA12878| 362581|    39922| 
| CHH_OT_NA12878| 420909|    44100| 
| CpG_CTOB_NA12878| 2407045|   721250| 
| CpG_CTOT_NA12878| 2420063|   728718| 
| CpG_OB_NA12878| 2901|     2200| 
| CpG_OT_NA12878| 3310|     2531|
