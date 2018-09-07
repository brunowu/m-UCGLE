#!/bin/bash

echo "=> Start cleaning ..."

rm -rf ${1%.*}"_cleaned"

mkdir ${1%.*}"_cleaned"

echo "=> ALL the cleaned results of $1 are seperately saved in the directory ${1%.*}"_cleaned""

grep "residual \[" $1 > tmp.txt

echo "=> Loading file: "${1}" ..."

cleanFile=${1%.*}"_clean.txt"

awk '{print $3 " " $6 }' tmp.txt > ${1%.*}"_cleaned"/${cleanFile}

BSIZE=`awk 'BEGIN { max=0 } $1 > max { max=$1;} END { print max }' FS=" " ${1%.*}"_cleaned"/${cleanFile}`

echo "=> There are "${BSIZE}" right hand sides in the same time in this file ..."

for i in `seq 0 ${BSIZE}`;
do
    echo "=>The cleaned results of RHS No. $i is saved into ${1%.*}"_cleaned"/${cleanFile}_RHS_$i.txt"
    grep "$i " ${1%.*}"_cleaned"/${cleanFile} | awk '{print $2}' > ${1%.*}"_cleaned"/${cleanFile}_RHS_$i.txt
done 

echo "=> Completed!"

rm tmp.txt ${1%.*}"_cleaned"/${cleanFile}
