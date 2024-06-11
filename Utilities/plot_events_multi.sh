#!/usr/bin/env bash

cd plot_events

f1=5_2_-1_1_2_-2_21_1_2_3_5_4
f2=5_2_-1_1_2_-2_21_1_4_3_5_2
f3=5_2_-1_1_2_-2_21_1_5_2_3_4
f4=5_2_-1_1_2_-2_21_1_5_4_3_2
#f5=6_2_2_3_0_7
#f6=6_2_2_2_1_7
#f7=6_2_2_1_2_7
#f8=6_2_2_0_3_7
#f9=6_2_2_2_0_0
#f10=6_2_2_0_1_1
#f11=6_2_2_1_0_1
#f12=6_2_2_1_1_0

bash plot_events.sh $f1
bash plot_events.sh $f2
bash plot_events.sh $f3
bash plot_events.sh $f4
#bash plot_events.sh $f5
#bash plot_events.sh $f6
#bash plot_events.sh $f7
#bash plot_events.sh $f8
#bash plot_events.sh $f9
#bash plot_events.sh $f10
#bash plot_events.sh $f11
#bash plot_events.sh $f12


##  2
cd ..
#./plot_events/internal/histograms.py Outputs/events__$f1.HwU Outputs/events__$f2.HwU --out=comparison --no_stat --assign_types='1-7-4-7','2-7-4-7' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact

##  6
#./plot_events/internal/histograms.py Outputs/events__$f1.HwU Outputs/events__$f2.HwU Outputs/events__$f3.HwU Outputs/events__$f4.HwU Outputs/events__$f5.HwU Outputs/events__$f6.HwU --out=comparison --no_stat --assign_types='1-0-0-1','1-0-1-0','1-1-0-0','2-0-0-1','2-0-1-0','2-1-0-0' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact

##  8

#./plot_events/internal/histograms.py Outputs/events__$f1.HwU Outputs/events__$f2.HwU Outputs/events__$f3.HwU Outputs/events__$f4.HwU Outputs/events__$f5.HwU Outputs/events__$f6.HwU Outputs/events__$f7.HwU Outputs/events__$f8.HwU --out=comparison --no_stat --assign_types='1-7-3-0','1-7-2-1','1-7-1-2','1-7-0-3','2-3-0-7','2-2-1-7','2-1-2-7','2-0-3-7' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact

### #12

#./plot_events/internal/histograms.py Outputs/events__$f1.HwU Outputs/events__$f2.HwU Outputs/events__$f3.HwU Outputs/events__$f4.HwU Outputs/events__$f5.HwU Outputs/events__$f6.HwU Outputs/events__$f7.HwU Outputs/events__$f8.HwU Outputs/events__$f9.HwU Outputs/events__$f10.HwU Outputs/events__$f11.HwU Outputs/events__$f12.HwU --out=comparison --no_stat --assign_types='1-0-0-2','1-0-2-0','1-2-0-0','1-0-1-1','1-1-0-1','1-1-1-0','2-0-0-2','2-0-2-0','2-2-0-0','2-0-1-1','2-1-0-1','2-1-1-0' --n_ratios=0 --no_scale --no_merging --no_pdf --no_alpsfact

#mv events_* Outputs

