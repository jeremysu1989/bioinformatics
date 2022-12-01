#!/bin/bash

function revcom()
{
 seq=$1
 revcom=`echo $seq| tr [AGCTagctUu] [TCGAtcgaAa] | rev`
 echo $revcom
}    
       
inputfile=$1

while read line
do
     seq="`echo $line | grep '^[AGCT]'`"
     if [ x"$seq" != x"" ]
     then
         echo $seq >> extracted.fasta
         revcom $seq >> revcom.extracted.fasta
     fi
done < $inputfile
