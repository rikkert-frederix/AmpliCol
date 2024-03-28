#! /usr/bin/env sh


gfortran -fcheck=all -o plot_events HwU.f analysis.f plot_events.f90

#rm -f events_*.HwU > /dev/null

i=0
for file in ../Outputs/events__$1_.lhe.rwgt ; do
    ./plot_events "$file"
    mv events.HwU events__$1.HwU
    ((i = i + 1))
done

#./internal/histograms.py events_*.HwU --average --no_stat  --out='../events__'$1 --no_open

#./internal/histograms.py events_*.HwU --sum --no_stat  --out='sum_2qq_5__'$1 --no_open


#mv ../events_* ../Outputs
