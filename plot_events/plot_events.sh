#! /usr/bin/env sh

gfortran -o plot_events plot_events.f90 analysis.f HwU.f
rm -f events_*.HwU > /dev/null

i=0
for file in ../events_*$1_2_?.lhe.rwgt ; do
    echo $file | ./plot_events
    mv events.HwU events_$i.HwU
    ((i = i + 1))
done

./internal/histograms.py events_*.HwU --sum --no_stat  --out='../events_'$1 --no_open
