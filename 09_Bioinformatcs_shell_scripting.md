a large part of bioinformatics is patching together various processing steps into a pipeline, and then repeatedly applying this pipeline to many files. this isn't exciting scientific work, but it's a necessary hurdle before tackling more exciting analyses.

how to wirte rerunnable bash shell scripts, automate file-processing tasks with *find* and *xarg*, run pipelines in parallel, and see a simple makefile

### Basic Bash Scripting
python may be a more suitable language for commonly reused or advanced pipelines

#### Variables and command arguments
Processing pipelines having numeroous settings that should be stored in variables (e.g., which directories to store results in, parameter values for commands, input files, etc.).  Bash also reads command-line arguments into variables, so you'll need to be familiar with accessing variables's values to work with command-line arguments.

To access a variables's value, we use a dollor sign in front of the variable's name. In some cases it's not clear where a variable name ends and where an adjacent string begins. To prevent this, wrap the variable name in braces. There is one more step we should take to make this more robust: quoting variables. If your script requires numerous or complicated options, it might be easier to use Python instead of Bash.
```
results_dir="results/"
echo "{$results_dir}"

sample="CONTROL"
mkdir "${sample}_aln"
```

#### Processing Files with Bash using for loops and globbing
in bioinformatics, most of our data is split across multiple files. At the heart of any processing pipeline is some way to apply the same workflow to each of these files, taking care to keep track of sample names.
```
#!/bin/bash
set -e
set -u
set -o pipefail

if [ "$#" -ne 1 ]
then
     echo "error:  you provided $#, 1 required"
     echo "usage: pipe_prep.sh samples.txt"
     exit 1
fi
# specify the input sample file, where the third column is the path to each sample FASTQ file
sample_info="$1"

# creat a Bash array from the third column of $sample_info
sample_file=($(cut -d ' ' -f 5 "$sample_info")) #subcommand substitution

mkdir stats

for fastq_file in ${sample_file[@]}
do
     # strip .fastq from each FASTQ file, and add suffix "-stats.txt" to create an output filename for each FASTQ file
     results_file="$(basename $fastq_file .fq)-stats.txt"
     # run fastq_stat on a file, writing results to the filename we've above
     wc -l $fastq_file > stats/$results_file
done
```
however, many bioinformatics pipelines combine two or more input fiels into a single output file. When writing script to algin paired-end reads, we can't loop over each file like we did earlier. Instead, each **sample**, rather than each **file**, is the processing unit. Consequently, our loop must iterate over unique samples names, and we use these samples names to re-create the input FASTQ fiels used in alignment.

```
#!/bin/bash
set -e
set -u
set -o pipefail

if [ "$#" -ne 1 ] 
then 
     echo "error:  you provided $#, 1 required"
     echo "usage: pipe_prep.sh samples.txt"
     exit 1
fi
sample_info="$1"

# get the file path
sample_path=$(dirname `head -n1 samples.txt | cut -d ' ' -f5`)

# get the sample name
sample_name=$((cut -d ' ' -f 1 "$sample_info") | uniq)

mkdir result

for sample in ${sample_name[@]}
do
     paste -d * $sample_path/${sample}_R0.fq $sample_path/${sample}_R1.fq > result/${sample}.res
done
```
Bash's loops are a handy way of applying commands to numerous files, but have a few downsides. First, compared to the Unix tool *find*, globbing is not a very powerful way to select certain files. Second, Bash's loop syntax is lengthy for simple operations, and a bit archaic. Finally, there's no easy way to parallelize Bash loops in a way that constrains the number of subprocesses used. 

#### Automating File-processing with *find* and *xargs*
using find and xargs
*find* 命令是Linux中强大的搜索命令，不仅可以按照文件名搜索文件，还可以按照权限文件大小、时间属性、inode等来搜索文件。

命令格式 
*find*命令的基本格式为： find 搜索路径 [选项] 搜索内容
*find*命令具有两个参数：
1. 第一个参数指定搜索路径
2. 第二个参数指定搜索内容

按照文件名搜索

支持的选项
1. -name 按照文件名搜索
2. -iname 按照文件名搜索，不区分文件名大小写
3. -inum 按照inode号搜索
```
find shell -name *.sh
find shell -iname *.SH
```

按照修改时间搜索

支持的选项
1. -atime [+-] 按照文件访问时间搜索
2. -mtine [+-] 按照文件数据修改时间搜索
3. -ctime [+-] 按照文件状态修改时间搜索

这里对+-的含义说明
1. -5：表示 5 天内修改的文件
2. 5：表示当前时间之前第 5-6 那一天修改的文件
3. +5：表示 6 天之前修改的文件
