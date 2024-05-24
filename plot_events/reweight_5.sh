#./plot_events.sh 5_2_-1_1_2_-2_21_1_2_3_5_4
#./plot_events.sh 5_2_-1_1_2_-2_21_1_5_2_3_4

./plot_events.sh 5_2_-1_1_2_-2_21_1_5_4_3_2
./plot_events.sh 5_2_-1_1_2_-2_21_1_4_3_5_2

./internal/histograms.py events_*.HwU --sum --no_stat  --out='sum_2qq_5__'$1 --no_open

