#!/bin/sh

# Needs #segs

clang -Xpreprocessor -fopenmp -lomp -w -o SpecFit SpecFit.c Utilities.c -lgsl  -lm

segs=$1

echo "# Segments =" $segs

for ((j=0;j<segs;j++))
do
echo $j
foo=$(printf "AET_seg%d_t.dat" $j)
./SpecFit $foo 0
foo=$(printf "specfit_0_%d.dat" $j)
mv specfit.dat $foo
gnuplot Qscan.gnu
foo=$(printf "Qscan_0_%d.png" $j)
mv Qscan.png $foo
foo=$(printf "AET_seg%d_t.dat" $j)
./SpecFit $foo 1
foo=$(printf "specfit_1_%d.dat" $j)
mv specfit.dat $foo
gnuplot Qscan.gnu
foo=$(printf "Qscan_1_%d.png" $j)
mv Qscan.png $foo
done


