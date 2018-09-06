#!/bin/bash

grep "residual " $1 > tmp.txt

echo "=> Loading file: "${1}" ..."

cleanFile=${1%.*}"_clean.txt"

awk '{print $3 " " $6 }' tmp.txt > tmp2.txt

BSIZE=`awk 'BEGIN { max=0 } $1 > max { max=$1;} END { print max }' FS=" " tmp2.txt`

echo "=> There are "${BSIZE}" right hand sides in the same time in this file ..."

awk 'BEGIN{print '$BSIZE'}{print}' tmp2.txt > ${cleanFile}

echo "!!> Attention, the first line of "${cleanFile}" gives the right hand sides number in the given file ..."
echo "=> The cleaned data is saved into : "${cleanFile}

echo "=> Completed!"
rm tmp.txt tmp2.txt


