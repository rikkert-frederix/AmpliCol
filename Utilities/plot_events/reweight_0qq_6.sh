#! /usr/bin/env sh


rm -f events_*.HwU > /dev/null

./plot_events.sh ../Outputs/0qq_randinit_40 6_2_21_21_21_21_21_21_1_2_3_4_5_6
./plot_events.sh ../Outputs/0qq_randinit_40 6_2_21_21_21_21_21_21_1_3_2_4_5_6
./plot_events.sh ../Outputs/0qq_randinit_40 6_2_21_21_21_21_21_21_1_3_4_2_5_6
./plot_events.sh ../Outputs/0qq_randinit_40 6_2_21_21_21_21_21_21_1_3_4_5_2_6
./plot_events.sh ../Outputs/0qq_randinit_40 6_2_21_21_21_21_21_21_1_3_4_5_6_2


./internal/histograms.py events_*.HwU --sum --no_stat  --out='sum_0qq_6__'$1 --no_open

