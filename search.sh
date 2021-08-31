#!/bin/sh

# Needs start and end segs

clang -Xpreprocessor -fopenmp -lomp -w -o search search.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas  -lm

segs=$1
sege=$2

echo "# start =" $segs
echo "# end =" $sege

for ((j=segs;j<=sege;j++))
do
echo $j
./search $j
done


