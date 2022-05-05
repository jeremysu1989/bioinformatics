# 将所有小写字母转换为大写字母
cat testfile | tr a-z A-Z
# -c 用字符串1中的字符集的补集替换字符串1； -s 删除重复的字符串； 用换行符替换转化后的第一个字符串，从而实现将有空格的文本分割成行
cat testfile | tr -cs 'a-zA-Z' '\n'
# 对上述命令进行组合即可统计文本文档中所有单词的出现频率,并输出前10个单词
cat testfile | tr -cs 'a-zA-Z' '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed 10q
# 同时查看文件的开头和结尾
(head -n2; tail -n2) < testfile
# 去除文件第一行，输出所有
tail -n +2 testfile
# les命令非常有用
less    
g 移动到首行
G 移动到末行
/ 向后匹配查询
? 向前匹配查询
