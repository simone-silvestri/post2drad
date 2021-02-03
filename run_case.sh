#!/bin/bash                                                                                                                                                                                                

echo $FILE0

numfiles=10

prefix=C 
cd /scratch-local/simone/POST_frag/data
dirlist=(`ls ${prefix}.*`)
length=${#dirlist[@]}
file=${dirlist[${length}-1]}
num=$(echo ${file} | sed 's/[^0-9]*//g')
num2=$(echo $num | sed 's/^0*//')
let "num2++"

if squeue -u simone | grep 'ppr.sh' &> logfile ; then
  echo "post still going"
else
   cd /scratch-local/simone/POST_frag/data/
   i=$num2
   echo $i
   startfile=$i
   if [[ $i -lt 330 ]]; then 
     rm U.* V.* W.* P.* C.*
     length=$(($numfiles+$i-1))
     echo $length
     while [[ $i -le $length ]]; do
         printf -v num2 "%05d" $i
         echo "copying number" $num2
         for j in C U V W P; do
             rsync -rv /archive/simone/backup-scratch/develop2/DATA/$j.$num2 .
         done
         ((i=i+1))
     done
     cd ../
     cp params_orig.f90 params.f90
     sed -i -e "s/SKIPFILES/$startfile/g" params.f90
     make clean 
     make 
     sbatch ppr.sh
     echo "restarted post"
   else
     echo "finished computations for average"
   fi
fi
echo "Finished testing"
