#!/bin/bash
# give the fastq file after the script name

inputfile=$1

while read line
do
     seq="`echo $line | grep '^[AGCT]'`"
     if [ x"$seq" != x"" ]
     then
         echo $seq >> extracted.fastq
     fi
done < $inputfile
