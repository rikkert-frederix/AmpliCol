#! /usr/bin/env sh

full_path='/home/timea/Documents/Uppsala_MUNKAHELY/Projects/Colour_implement/IntegrateGluons/Utilities/plot_events'

cd $full_path

gfortran -fcheck=all -o plot_events HwU.f analysis.f plot_events.f90
echo $2
mkdir "res_wgt_$2"
#rm -f events_*.HwU > /dev/null

for file in $1_*.lhe.rwgt ; do
	./plot_events $file
	mv events.HwU "res_wgt_$2"/events_$3.HwU
done

cd res_wgt_$2
../internal/histograms.py events_*.HwU --sum --no_stat  --out='sum_'$2 --no_open
