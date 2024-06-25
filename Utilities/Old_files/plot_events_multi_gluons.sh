#!/usr/bin/env bash

cd plot_events

bash plot_events_gluons.sh 4
bash plot_events_gluons.sh 5
#bash plot_events_gluons.sh 6
#bash plot_events_gluons.sh 7
#bash plot_events_gluons.sh 8

#cd ..

#./plot_events/internal/histograms.py Outputs/events_4_2_0_0_2_0.HwU Outputs/events_5_2_0_0_3_0.HwU Outputs/events_6_2_0_0_4_0.HwU Outputs/events_7_2_0_0_5_0.HwU Outputs/events_8_2_0_0_6_0.HwU --out=comparison --no_stat --assign_types='n=4' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact

#mv events_* Outputs

