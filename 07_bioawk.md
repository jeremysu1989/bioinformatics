## awk的升级版bioawk
在awk的基础上添加以下功能
- 将FASTA / FASTQ文件看作单行来处理
- 为已知数据格式添加了新的内部变量。
- 添加了许多与生物信息学相关的功能，例如revcomp，以产生反向互补序列

#### bioawk的安装
```
conda install bioawk
#or
brew install bioawk
```

#### bioawk的帮助文档
```
bioawk -c help
bed:
	1:chrom 2:start 3:end 4:name 5:score 6:strand 7:thickstart 8:thickend 9:rgb 10:blockcount 11:blocksizes 12:blockstarts 
sam:
	1:qname 2:flag 3:rname 4:pos 5:mapq 6:cigar 7:rnext 8:pnext 9:tlen 10:seq 11:qual 
vcf:
	1:chrom 2:pos 3:id 4:ref 5:alt 6:qual 7:filter 8:info 
gff:
	1:seqname 2:source 3:feature 4:start 5:end 6:score 7:filter 8:strand 9:group 10:attribute 
fastx:
	1:name 2:seq 3:qual 4:comment

man bioawk # for more help
```
