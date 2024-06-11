./plot_events.sh 6_2_-1_1_2_-2_21_21_1_2_3_5_6_4
./plot_events.sh 6_2_-1_1_2_-2_21_21_1_2_3_6_5_4
./plot_events.sh 6_2_-1_1_2_-2_21_21_1_5_2_3_6_4
./plot_events.sh 6_2_-1_1_2_-2_21_21_1_5_6_2_3_4
./plot_events.sh 6_2_-1_1_2_-2_21_21_1_6_2_3_5_4
./plot_events.sh 6_2_-1_1_2_-2_21_21_1_6_5_2_3_4

./plot_events.sh 6_2_-1_1_2_-2_21_21_1_4_3_5_6_2
./plot_events.sh 6_2_-1_1_2_-2_21_21_1_4_3_6_5_2
./plot_events.sh 6_2_-1_1_2_-2_21_21_1_5_4_3_6_2
./plot_events.sh 6_2_-1_1_2_-2_21_21_1_5_6_4_3_2
./plot_events.sh 6_2_-1_1_2_-2_21_21_1_6_4_3_5_2
./plot_events.sh 6_2_-1_1_2_-2_21_21_1_6_5_4_3_2




./internal/histograms.py events_*.HwU --sum --no_stat  --out='sum_2qq_6__'$1 --no_open

