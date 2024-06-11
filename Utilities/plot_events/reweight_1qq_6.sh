#! /usr/bin/env sh


rm -f events_*.HwU > /dev/null


./plot_events.sh ../Outputs/1qq_qq 6_2_-1_1_21_21_21_21_1_3_4_5_6_2



./internal/histograms.py events_*.HwU --sum --no_stat  --out='sum_1qq_6__'$1 --no_open

