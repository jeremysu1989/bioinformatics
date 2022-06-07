## 文件内容操作 
#### 将所有小写字母转换为大写字母
    cat testfile | tr a-z A-Z
#### -c 用字符串1中的字符集的补集替换字符串1； -s 删除重复的字符串； 用换行符替换转化后的第一个字符串，从而实现将有空格的文本分割成行
    cat testfile | tr -cs 'a-zA-Z' '\n'
#### 对上述命令进行组合即可统计文本文档中所有单词的出现频率,并输出前10个单词
    cat testfile | tr -cs 'a-zA-Z' '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed 10q
#### 同时查看文件的开头和结尾
    (head -n2; tail -n2) < testfile
#### 去除文件第一行，输出所有
    tail -n +2 testfile
#### less命令非常有用
less    
g 移动到首行
G 移动到末行
/ 向后匹配查询
? 向前匹配查询
#### 对于tab分隔的文件，选择性的展示某一列
    cut -f 2 testfile
#### 使用cut命令，结合-d参数指明delimiter character可以进行分割,例如使用大写S作为分割的字符
```
cut -dS -f3 testfile
    
# 读取fastq文件的ID序列，用：进行区域分隔，并输出第4-6和7区域的内容，现实前5行
bioawk -c fastx '{print $1}' NA12878_R0.fq | cut -d : -f4-6,7 | head -n5

# 将 PATH 变量取出，我要找出第二到第四，还有第六个路径，并分行展示
echo $PATH | cut -d ":" -f 2-4,6 | tr ":" "\n"
```
#### 对于空格分隔的文件，选择性的展示某一列
    awk -F " " '{print $2}' testfile
#### 使用colomn命令实现文本文件的表格化展示，-s参数可以设定分隔字符
    column -s "," -t testfile
#### grep命令对搜索字符的精确匹配，按照单词进行匹配
    grep -w "green" testfile
#### 可以通过添加正则表达的方式让搜索命令更强大
    grep -w -E "as|be|the" testfile
#### 添加-c可以统计实现匹配的行数
    grep -c -w -E "as|the" testfile
#### -o只输出匹配的字符
    grep -o -E "as|the" testfile
#### \w 匹配字母或数字或下划线或汉字 等价于 '[^A-Za-z0-9]',“+”元字符规定其前导字符必须在目标对象中连续出现一次或多次
    grep -E -o 'gene_id "(\w+)"' testfile
#### 根据上面的命令，管道连接几个命令可以提取gtf文件中所有基因的名称
    grep -E -o 'gene_id "(\w+)"' Sars_cov_2.ASM985889v3.101.gtf | cut -f2 -d' ' | sort | uniq | sed 's/"//g'
#### 对列表文件按照第一列，第二列进行排序sort
    sort -k1,1 -k2,2 testfile
#### 对列表文件按照第一列排序，第二列按照数字排序
    sort -k1,1 -k2,2n testfile
#### 对上述排序结果进行倒序输出
    sort -k1,1 -k2,2n -r testfile
#### 仅对排序后的第二列结果进行倒序输出
    sort -k1,1 -k2,2nr testfile
#### 经常使用的组合
    sort | uniq   or  sort | uniq -c
#### 查看gtf文件中基因组上各个基因原件的个数
    grep -v "^#" Sars_cov_2.gtf | cut -f3 | sort | uniq -c | sort -rn

Without having to load data into a program like R or Excel, we can quickly calculate
summary statistics about our plain-text data files.

---

#### awk进行文本文件的处理,模式匹配的条件下执行动作
    awk pattern { action } testfile
#### awk输出符合pattern的records
    awk pattern testfile
#### 输出位于chr1上面且碱基个数大于10的records
    awk "$1 ~ /chr1/ && $3-S2 > 10" testfile
#### awk直接对records进行操作
    awk {action} testfile
#### 输出位于chr2或chr3上的records，并且计算该records上碱基的数目
    awk "$1 ~ /chr2|chr3/" {print $0 "\t" $3-$2} testfile
#### 将csv文件转换成tab分隔文件
    awk -F"," -v OFS="\t" {print $1,$2,$3,$4} testfile.csv
#### 使用NR（number of field）选取特定行的内容
    awk "NR >=2 && NR <= 5" testfile
#### gtf文件转化为bed文件
    awk '!/^#/ { print $1 "\t" $4-1 "\t" $5 }' testfile
#### awk求某一列数字的和
    awk '{sum+=$5} END{print sum}'
#### awk计算指定行的和
    awk 'NR==3{for(i=1;i<=NF;i++)sum=sum+$i;}END{print sum}'
#### awk还有很多自带的函数功能，可以慢慢研究

#### echo输出结果时取消自动换行
```
echo -ne $i; # used in iteration
```
---

### sed在生物信息学数据处理中的常用命令
    sed 'action' testfile
    sed 's/pattern/replace/' testfile
#### 修改数据中chrom为chr，输出结果，并现实前3行
    sed 's/chrom/chr/g' testfile | head -n3
#### 忽略匹配字符的大小写
    sed 's/chrom/chr/i' testfile
#### sed对模式进行匹配操作，可以大大增强sed的功能
    echo "chr1:28427874-28425431" | sed -E 's/^(chr[^:]+):([0-9]+)-([0-9]+)/\1\t\2\t\3/'
##### —E 参数可以调用扩展的正则匹配方式；正则表达式的第一部分^(chr[^:]+):从行首开始匹配chr，匹配多个非冒号的字符，一直匹配到冒号；正则表达式的第二部分和第三部分([0-9]+)都是进行多个数字的匹配，第二三部分之间用‘-’连接；
#### 有其他更多的方式可以实现上述功能
    echo "chr1:28427874-28425431" | sed 's/[:-]/\t/g'
#### 或者
    echo "chr1:28427874-28425431" | sed 's/:/\t/' | sed 's/-/\t/'
#### 或者
    echo "chr1:28427874-28425431" | tr ':-' '\t'
#### 提取gtf文件中所有转录本的名称
    grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | head -n 3 | sed -E -n 's/.*transcript_id "([^"]+)".*/\1/p'
##### sed命令里面匹配的模式为transcript_id之前的任意多个任意字符，匹配transcript_id “，匹配大于等于一个非双引号字符，匹配”，匹配任意多个任意字符，然后输出匹配到的字符；-n参数使得sed不能输出任何行内容；-p参数确保只有匹配的行被输出
#### 使用sed命令输出文件的特定行
    sed -n "20,50p" testfile
#### subshell, 实现对gtf文件的排序
    (grep "^#" Sars_cov_2.gtf; grep -v "^#" Sars_cov_2.gtf | sort -k1,1 -k4,4n) > Sars_cov_2.sort.gtf



## 文件整体操作
#### 统计文件中行数
    wc -l
#### 统计文件中非空行数,
    grep -c "[^ //n//t]" testfile
#### 统计文件的列数
    awk -F "\t" '{print NF; exit}' testfile
#### 统计文件中符合条件的列数
    grep -v "^#" testfile | awk -F "\t" '{print NF; exit}'
#### 压缩文件的操作
 1. *.tar 用 tar –xvf解压 
 2. *.gz 用 gzip -d或者gunzip解压 
 3. *.tar.gz和*.tgz 用 tar –xzf解压 
 4. *.bz2 用 bzip2 -d或者用bunzip2解压 
 5. *.tar.bz2用tar –xjf解压 
 6. *.Z 用 uncompress解压 
 7. *.tar.Z 用tar –xZf解压 
 8. *.rar 用 unrar解压
 9. *.zip 用 unzip解压
#### 批量修改文件名称1 rename,类似于perl脚本的方式
    rename 's/abc/def/' *
###### 使用rename命令之前在macOS安装该命令
    brew update-reset
    brew --version
    brew install rename
#### wget下载单个文件至指定文件夹
    wget -p /Users/scarecrow/Temp/test/download -c website/filename
#### wget批量下载多个文件至指定文件夹，并切-b让程序在后台运行,这个命令提高下载速度（大约不到5s），如果没有-b参数，则wget会逐个文件进行解析，逐个下载，下载时间大约20分钟
    for i in $(awk '{print $3}' testfile); do wget -bc -P /Users/scarecrow/Temp/test/test_download http://ftp.ensembl.org/pub/release-106/fasta/saccharomyces_cerevisiae/dna/$i; done
#### wget将下载信息写入文件,将下载记录写入文件，这种会让下载界面很清爽，没有特别多的wget-log文件,耗时3s
    for i in $(awk '{print $3}' test.txt); do wget -o download.log -bc -P /Users/scarecrow/Temp/test/test_download http://ftp.ensembl.org/pub/release-106/fasta/saccharomyces_cerevisiae/dna/$i; done
    for i in $(awk '{print $3}' CHECKSUMS); do wget -o download.log -bc -P /Users/scarecrow/Temp/test/test_download http://ftp.ensembl.org/pub/release-106/fasta/saccharomyces_cerevisiae/dna/$i; done 
#### find命令及复制操作
    find . -name "*QC.txt" -exec cp '{}' /target/folder/ \;
#### 打印文件夹下文件的全路径
```
ls -lrt -d -1 "$PWD"/{*,.*}
ls | sed "s:^:$PWD/:"
ls | sed "s:^:`pwd`/:"
realpath filename
```


## 文件间操作
#### join函数，将两个不同的文件按照相同的元素进行拼接，拼接之前必须进行分别排序
    sort -k1,1 testfile1
    sort -k1,1 testfile2
    join -1 1 -2 1 testfile1 testfile2
#### 有名管道（FIFO）和无名管道（pipe）
https://blog.csdn.net/qq_41945905/article/details/110489607?utm_medium=distribute.pc_relevant.none-task-blog-2~default~baidujs_baidulandingword~default-4.pc_relevant_antiscanv2&spm=1001.2101.3001.4242.3&utm_relevant_index=7

https://blog.csdn.net/yxtxiaotian/article/details/69568774?spm=1001.2101.3001.6650.4&utm_medium=distribute.pc_relevant.none-task-blog-2%7Edefault%7ECTRLIST%7ERate-4.pc_relevant_antiscanv2&depth_1-utm_source=distribute.pc_relevant.none-task-blog-2%7Edefault%7ECTRLIST%7ERate-4.pc_relevant_antiscanv2&utm_relevant_index=9

#### paste函数，两个文件按行合并，按列合并，必须确保两个文件的行数相等
paste -d* file1 file2 可以让两个文件合并的时候中间插入*

paste -s file1 file2 可以让两个文件合并时生成两行
```
wc -l * | sed '$d' > file_line.txt
awk '{print $1 $2}' CHECKSUM > checksum.txt
wc -l file_line.txt checksum.txt
# 55 file_line.txt
# 55 checksum.txt
# 110 total
paste file_line.txt checksum.txt
     # 121 README	117755
     # 335 Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa.gz	3326870
     # 1075 Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.II.fa.gz	13070249
     # 403 Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.III.fa.gz	4414597
     # ...
     # total 55

paste file_line.txt checksum.txt > result.txt
awk '{print "\|"$3"\|"$2"\|"$1"\|"}' result.txt
```
You can use this result to generate table
|Checksum |File_name|Line_No|
|-----|-----|-----|
|117755|README|121|
|3326870|Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa.gz|335|
|13070249|Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.II.fa.gz|1075|
|4414597|Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.III.fa.gz|403|
|05312465|Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.IV.fa.gz|2008|




## 服务器和本地之间文件传输
所有的请求都是从本地电脑发出的
#### download file over SSH protocol
    scp username@hostname:/path/to/remote/file /path/to/local/file
    scp jsu@192.168.202.180:/data/home/jeremysu/test/fastq_checksums.sha /Users/scarecrow

#### Uploading a file from a local computer to a remote one
    scp /path/to/local/file username@hostname:/path/to/remote/file
    scp /Users/scarecrow/Temp/test/Sars_cov_2.gff3 jsu@192.168.202.180:/data/home/jeremysu/test/


## 终端进程管理
#### 批量终止user下的进程
```
 1609  ps x
 1610  ps x|cut -f 1 -d " "
 1611  ps x|cut -f 21 -d " "
 1612  ps x|cut -f 2 -d " "
 1613  kill `ps x|cut -f 2 -d " "|tr "\n" " "`
 ```
#### 终止用户名下所有进程
```
pkill -u jeremysu
killall -u jeremysu
```
#### 依次终止用户名下所有进程
```
ps -ef | grep "jeremysu" | awk '{print $2}' | sudo xargs kill -9
pgrep -u "jeremysu" | sudo xargs kill -9
```
