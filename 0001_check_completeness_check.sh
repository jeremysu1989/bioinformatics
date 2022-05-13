```
#!/bin/bash
file="$1"
dir="$2"

if [ -e $file ]
then
     echo "$file is exist!"
     if [ -d $dir ]
     then 
          echo "$dir is directory and exist!"
          if [ -f "download.log" ]
          then
               rm download.log
          fi  
          ######
          for i in $(ls $dir)
          do  
               a=`grep $i $file | awk '{print $1 $2}' | sed 's/^0//'`
               b=`(sum $dir/$i | awk '{print $1 $2}')`
               if [ $a = $b ]
               then
                    echo "$i is downloaded completelt, Congratulation!" >> download.log
               else
                    echo "Something wrong with $i, pls check manually!" >> download.log
               fi  
          done
          date >> download.log
          echo "Check finished!"
     else 
          echo "$dir is not exist!"
          break
     fi  
else
     echo "$file is not exist! Quit!"
     break
fi
date
```
