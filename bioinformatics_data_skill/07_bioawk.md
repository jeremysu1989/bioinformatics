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
### 处理sam文件
#### 查看sam文件中每个碱基的甲基化状态
```
bioawk -c sam '{print "seq :"$seq "\n" $14}' NA12878.sam | head -n8

seq :ACAAATCTATCACCCTATTAACCACTCACAAAAACTCTCCATACATTTAATATTTTCATCTAAAAAATATACACACAATAACATTACAAAACACTAAAACCAAAACACCC
XM:Z:...xh........................zxh.h........h.....hh.......z...xhhhhh.h.h...z.z...h....h.z.h..z..xh.h..zx.h.....
seq :TCACCCTATTAACCACTCACAAAAACTCTCCATACATTTAATATTTTCATCTAAAAAATATACACACAATAACATTACAAAACACTAAAACCAAAACACCACTA
XM:Z:....................zxh.h........h.....hh.......z...xhhhhh.h.h...z.z...h....h.z.h..z..xh.h..zx.h........
seq :CACCCTATTAACCACTCACAAAAACTCTCCATACATTTAATATTTTCATCTAAAAAATATACACACAATAACATTACAAAACACTAAAACCAAAACACCCTATATCACA
XM:Z:...................zxh.h........h.....hh.......z...xhhhhh.h.h...z.z...h....h.z.h..z..xh.h..zx.h........h..z..
seq :ACCCTATTAACCACTCACGAAAACTCTCCATACATTTAATATTTTCATCTAGAAAATATACACACAATAACATTACAAAACACTAAAACCAAAACACCC
XM:Z:..................Zxh.h........h.....hh.......z...xHhhhh.h.h...z.z...h....h.z.h..z..xh.h..zx.h.....
```

#### 使用if判断，选取序列为正向比对的reads,可以准确查看该位点的甲基化情况
```
bioawk -c sam '{if($2==0) print "seq :"$seq "\n" $14}' NA12878.sam | head -n8

seq :GGTGGGAGTTTTTTATGTATTTGGTATTTTTGTTTGGGGGGTATGTATGTGATAGTATTGTGAGATGTTGGAGTTGGAGTATTTTATGTTGTAGTATTTGTTTTTGATTTTTGTTTTATTTTATTATTTATTGTATTTATGTTTAATATTAT
XM:Z:..z.....h.h.hh...h............z..x...........h.z.z.....h....z....z.x.....xz....h.hhh.....z.x.....x...h.......hx..hh.h...h..........z.h.hh..z...h.......x
seq :GGGGTTGTATTTTTGTTTGGGGGGTATGTATGTGATAGTATTGTGAGATGTTGGAGTTGGAGTATTTTATGTTGTAGTATTTGTTTTTGATTTTTGTTTTATTTTATTATTTATTGTATTTATGTTTAATATTATAGGTGAATATATTTATTAAAGTGTGTTAATTAATTAATGT
XM:Z:.............z..x...........h.z.z.....h....z....z.x.....xz....h.hhh.....z.x.....x...h.......hx..hh.h...h..........z.h.hh..z...h.......x...z...h...hh..h.......................h
seq :GGGGTTGTATTTTTGTTTGGGGGGTATGTATGTGATAGTATTGTGAGATGTTGGAGTTGGAGTATTTTATGTTGTAGTATTTGTTTTTGATTTTTGTTTTATTTTATTATTTATTGTATTTATGTTTAATATTATAGGTGAATATATTTATTAAAGTGTGTTAATTAATTAATGTTTGTAGGATATAATAATAATAATTGAATGTTTGTATAGT
XM:Z:.............z..x...........h.z.z.....h....z....z.x.....xz....h.hhh.....z.x.....x...h.......hx..hh.h...h..........z.h.hh..z...h.......x...z...h...hh..h.......................h........h..........h..........x..h.x..x
seq :TTTGGTATTTTTGTTTGGGGGGTATGTATGTGATAGTATTGTGAGATGTTGGAGTTGGAGTATTTTATGTTGTAGTATTTGTTTTTGATTTTTGTTTTATTTTATTATTTATTGTATTTATGTTTAATATTATAGGTGAATATATTTATTAAAGTGTGTTAATTAATTAATGTTTG
XM:Z:...........z..x...........h.z.z.....h....z....z.x.....xz....h.hhh.....z.x.....x...h.......hx..hh.h...h..........z.h.hh..z...h.......x...z...h...hh..h.......................h...
```

#### 添加其他field数据的输出,使用length()输出序列长度
```
bioawk -c sam '{if($2==0) print length($seq)"bp" "\n" "seq :" $seq "\n" $14}' NA12878.sam | head -n9

152bp
seq :GGTGGGAGTTTTTTATGTATTTGGTATTTTTGTTTGGGGGGTATGTATGTGATAGTATTGTGAGATGTTGGAGTTGGAGTATTTTATGTTGTAGTATTTGTTTTTGATTTTTGTTTTATTTTATTATTTATTGTATTTATGTTTAATATTAT
XM:Z:..z.....h.h.hh...h............z..x...........h.z.z.....h....z....z.x.....xz....h.hhh.....z.x.....x...h.......hx..hh.h...h..........z.h.hh..z...h.......x
175bp
seq :GGGGTTGTATTTTTGTTTGGGGGGTATGTATGTGATAGTATTGTGAGATGTTGGAGTTGGAGTATTTTATGTTGTAGTATTTGTTTTTGATTTTTGTTTTATTTTATTATTTATTGTATTTATGTTTAATATTATAGGTGAATATATTTATTAAAGTGTGTTAATTAATTAATGT
XM:Z:.............z..x...........h.z.z.....h....z....z.x.....xz....h.hhh.....z.x.....x...h.......hx..hh.h...h..........z.h.hh..z...h.......x...z...h...hh..h.......................h
214bp
seq :GGGGTTGTATTTTTGTTTGGGGGGTATGTATGTGATAGTATTGTGAGATGTTGGAGTTGGAGTATTTTATGTTGTAGTATTTGTTTTTGATTTTTGTTTTATTTTATTATTTATTGTATTTATGTTTAATATTATAGGTGAATATATTTATTAAAGTGTGTTAATTAATTAATGTTTGTAGGATATAATAATAATAATTGAATGTTTGTATAGT
XM:Z:.............z..x...........h.z.z.....h....z....z.x.....xz....h.hhh.....z.x.....x...h.......hx..hh.h...h..........z.h.hh..z...h.......x...z...h...hh..h.......................h........h..........h..........x..h.x..x
```
#### 输出满足多个条件的sam序列，并转化成fasta
```
bioawk -c sam '{if($flag==16 && length($seq)<50) print ">" $qname "\n" revcomp($seq)}' NA12878.sam | head -n10
>A00265:368:H3M7FDSXY:3:2457:20491:34194:GATATCCATACTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT
TGTGATATAGGGTGTTTTGGTTTTAGTGTTTTGTAATGTTAT
>A00265:368:H3M7FDSXY:3:2541:13313:5791:ACAAGTCGGTGGCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT
TGATGTTTGTGTGGAAAGTGGTTGTGTAGATATTTAATTG
>A00265:368:H3M7FDSXY:3:1555:26648:18020:AATATTGAGGACCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT
GTGGTTAGAAGTGGGGGGAGGGGGGGGTTTGGTGGAAATTTTTTGTTA
>A00265:368:H3M7FDSXY:3:1645:31955:4883:TGATGTACATACCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT
GTGGTTAGAAGTGGGGGGAGGGGGGGGTTTGGTGGAAATTTTTTGTT
>A00265:368:H3M7FDSXY:3:1677:28528:21872:AAGGTCTTTAAACTTCTTAT_1:N:0:ACACTAAG+ATCCATAT
TTTGGGGTTTGGTAGAGATGTGTTTAAGTGTTGTGGTTAGAAGTGGGGG
```
#### Get the mean Phred quality score
```
bioawk -c sam '{print $qname, meanqual($qual)}' NA12878.sam | head -n5
A00265:368:H3M7FDSXY:3:1448:32280:17910:GAGTAGTTCCGTCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT	39.7545
A00265:368:H3M7FDSXY:3:2613:19578:6934:GGGTCGTTTGTCCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT	39.7885
A00265:368:H3M7FDSXY:3:2534:15492:31078:AGGTCCCGCAAACTTCTTAT_1:N:0:ACACTAAG+ATCCATAT	40
A00265:368:H3M7FDSXY:3:2451:9570:31876:ACCGTGAAATCACTTCTTAT_1:N:0:ACACTAAG+ATCCATAT	39.9697
A00265:368:H3M7FDSXY:3:2415:12237:5494:CCAAATTCTCGCCTTCTTAT_1:N:0:ACACTAAG+ATCCATAT	39.5263
```

### 处理fasta文件
#### 截取每一段序列的长度
```
bioawk -c fastx '{print $name, length($seq), $3}' Sars_cov_2.dna.toplevel.fa 

MN908947.3	29903
```

#### 测量fasta序列文件中gc含量
```
bioawk -c fastx '{print $name, length($seq), gc($seq)}' Sars_cov_2.dna.toplevel.fa 
MN908947.3	29903	0.379728
```

#### 获取所有序列的反向互补链序列
```
bioawk -c fastx '{print $seq, revcomp($seq)}' Sars_cov_2.dna.toplevel.fa 
```

#### convert fasta to tabular format
```
bioawk -t -c fastx '{ print $name, $seq }' input.fasta
```

### For fastq
Here, the -c fastx option remains same but bioawk will automatically recognize the fastq format and build the required variables, such as $name $seq $qual and $comment
```
#查看fastq文件中每条reads碱基的平均测序质量
bioawk -c fastx '{print meanqual($qual)}' NA12878_R0.fq | head -n10

#展示fastq文件中reads的ID
bioawk -c fastx '{print $name}' NA12878_R0.fq | head -n5
```
### for bed files
```
bioawk -c bed '{ print $end - $start }' test.bed
```
