#!/usr/bin/env bash

cd plot_events
./plot_events.sh 4
./plot_events.sh 5
./plot_events.sh 6
./plot_events.sh 7
./plot_events.sh 8
./plot_events.sh 9
#./plot_events.sh 10
#./plot_events.sh 11
#./plot_events.sh 12
cd ..
#./plot_events/internal/histograms.py  events_4.HwU events_5.HwU events_6.HwU events_7.HwU events_8.HwU events_9.HwU events_10.HwU events_11.HwU events_12.HwU --out=comparison --no_stat --assign_types='n=4','n=5','n=6','n=7','n=8','n=9','n=10','n=11','n=12' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact
#./plot_events/internal/histograms.py  events_4.HwU events_5.HwU events_6.HwU events_7.HwU events_8.HwU events_9.HwU events_10.HwU events_11.HwU --out=comparison --no_stat --assign_types='n=4','n=5','n=6','n=7','n=8','n=9','n=10','n=11' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact
#./plot_events/internal/histograms.py  events_4.HwU events_5.HwU events_6.HwU events_7.HwU events_8.HwU events_9.HwU events_10.HwU --out=comparison --no_stat --assign_types='n=4','n=5','n=6','n=7','n=8','n=9','n=10' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact
./plot_events/internal/histograms.py events_4.HwU events_5.HwU events_6.HwU events_7.HwU events_8.HwU events_9.HwU --out=comparison --no_stat --assign_types='n=4','n=5','n=6','n=7','n=8','n=9' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact
#./plot_events/internal/histograms.py events_4.HwU events_5.HwU events_6.HwU events_7.HwU events_8.HwU --out=comparison --no_stat --assign_types='n=4','n=5','n=6','n=7','n=8' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact
#./plot_events/internal/histograms.py events_4.HwU events_5.HwU events_6.HwU events_7.HwU --out=comparison --no_stat --assign_types='n=4','n=5','n=6','n=7' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact
#./plot_events/internal/histograms.py events_4.HwU events_5.HwU events_6.HwU --out=comparison --no_stat --assign_types='n=4','n=5','n=6' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact
