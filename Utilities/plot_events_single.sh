#!/usr/bin/env bash

cd plot_events

bash ./plot_events_single.sh 5_2_0_0_0_0


cd ..
#./plot_events/internal/histograms.py Outputs/events_4_2_0.HwU --out=comparison --no_stat --assign_types='n=4' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact


mv events_* Outputs

