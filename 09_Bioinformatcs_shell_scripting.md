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
