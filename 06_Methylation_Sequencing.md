## Analysing bisulfite methylation sequencing data

---

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
