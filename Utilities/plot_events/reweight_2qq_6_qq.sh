#! /usr/bin/env sh

all=1

folder=../Outputs/2qq_qq

if [ $all -eq 1 ] 
then

rm -f events_*.HwU > /dev/null
./plot_events.sh $folder 6_2_-1_1_2_-2_21_21_1_4_3_5_6_2
./plot_events.sh $folder 6_2_-1_1_2_-2_21_21_1_5_4_3_6_2
./plot_events.sh $folder 6_2_-1_1_2_-2_21_21_1_5_6_4_3_2

./plot_events.sh $folder 6_2_-1_1_2_-2_21_21_1_2_3_5_6_4
./plot_events.sh $folder 6_2_-1_1_2_-2_21_21_1_5_2_3_6_4
./plot_events.sh $folder 6_2_-1_1_2_-2_21_21_1_5_6_2_3_4
./internal/histograms.py events_*.HwU --sum --no_stat  --out='sum_2qq_6_qq_all_'$1 --no_open
fi

if [ $all -eq 2 ]
then

rm -f events_*.HwU > /dev/null
./plot_events.sh $folder 6_2_-1_1_2_-2_21_21_1_4_3_5_6_2
./plot_events.sh $folder 6_2_-1_1_2_-2_21_21_1_5_4_3_6_2
./plot_events.sh $folder 6_2_-1_1_2_-2_21_21_1_5_6_4_3_2
./internal/histograms.py events_*.HwU --sum --no_stat  --out='sum_2qq_6_qq_LC_'$1 --no_open
fi
