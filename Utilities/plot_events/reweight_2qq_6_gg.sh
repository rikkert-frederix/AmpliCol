#! /usr/bin/env sh

all=1

folder=../Outputs/2qq_gg_100k

if [ $all -eq 1 ] 
then

rm -f events_*.HwU > /dev/null
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_1_2_4_5_6
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_2_1_4_5_6
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_1_4_5_2_6
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_2_4_5_1_6
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_4_5_2_1_6
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_4_5_2_1_6

./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_1_2_6_5_4
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_2_1_6_5_4
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_1_6_5_2_4
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_2_6_5_1_4
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_6_5_1_2_4
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_6_5_2_1_4

./internal/histograms.py events_*.HwU --sum --no_stat  --out='sum_2qq_6_gg_all_'$1 --no_open
fi

if [ $all -eq 2 ]
then

rm -f events_*.HwU > /dev/null
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_1_2_6_5_4
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_2_1_6_5_4
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_1_6_5_2_4
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_2_6_5_1_4
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_6_5_1_2_4
./plot_events.sh $folder 6_2_21_21_1_-1_2_-2_3_6_5_2_1_4
./internal/histograms.py events_*.HwU --sum --no_stat  --out='sum_2qq_6_gg_LC_'$1 --no_open
fi
