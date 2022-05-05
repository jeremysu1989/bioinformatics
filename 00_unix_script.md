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
cut -dS -f3 testfile
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
 1、*.tar 用 tar –xvf解压 
 2、*.gz 用 gzip -d或者gunzip解压 
 3、*.tar.gz和*.tgz 用 tar –xzf解压 
 4、*.bz2 用 bzip2 -d或者用bunzip2解压 
 5、*.tar.bz2用tar –xjf解压 
 6、*.Z 用 uncompress解压 
 7、*.tar.Z 用tar –xZf解压 
 8、*.rar 用 unrar解压
 9、*.zip 用 unzip解压
#### 批量修改文件名称1 rename,类似于perl脚本的方式
rename 's/abc/def/' *
###### 使用rename命令之前在macOS安装该命令
brew update-reset
brew --version
brew install rename
m

