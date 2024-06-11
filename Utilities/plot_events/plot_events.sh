#! /usr/bin/env sh


gfortran -fcheck=all -o plot_events HwU.f analysis.f plot_events.f90

#rm -f events_*.HwU > /dev/null

full_path='/home/timea/Documents/Uppsala_MUNKAHELY/Projects/Colour_implement/IntegrateGluons-2qq_gluons/plot_events'
i=0
cd $1
for file in events__$2_.lhe.rwgt ; do
    $full_path/plot_events "$file"
    mv events.HwU $full_path/events__$2.HwU
    ((i = i + 1))
done

#./internal/histograms.py events_*.HwU --average --no_stat  --out='../events__'$1 --no_open

#./internal/histograms.py events_*.HwU --sum --no_stat  --out='sum_2qq_5__'$1 --no_open

#./internal/histograms.py sum_0qq_6__.HwU sum_1qq_6__.HwU sum_2qq_6_all_.HwU sum_2qq_6_LC_.HwU --out=comparison --no_stat --assign_types='0qq','1qq','2qq ALL','2qq LC' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact

#./internal/histograms.py sum_0qq_6__.HwU --out=kinem_0qq --no_stat
#./internal/histograms.py sum_1qq_6__.HwU --out=kinem_1qq --no_stat
#./internal/histograms.py sum_2qq_6_all_.HwU --out=kinem_2qq_all --no_stat
#./internal/histograms.py sum_2qq_6_LC_.HwU --out=kinem_2qq_LC --no_stat


#mv ../events_* ../Outputs
