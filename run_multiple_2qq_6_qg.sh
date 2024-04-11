#! /usr/bin/env sh

make 
# 2

./matrix_integrate_QCD 1 0 6 -1 21 -1 2 -2 21 1 2 6 3 4 5  > output_1_0_6_-1_21_-1_2_-2_1_2_6_3_4_5.txt
./matrix_integrate_QCD 1 1 6 -1 21 -1 2 -2 21 1 2 6 3 4 5  > output_1_1_6_-1_21_-1_2_-2_1_2_6_3_4_5.txt
./matrix_integrate_QCD 1 2 6 -1 21 -1 2 -2 21 1 2 6 3 4 5  > output_1_2_6_-1_21_-1_2_-2_1_2_6_3_4_5.txt

./matrix_integrate_QCD 1 0 6 -1 21 -1 2 -2 21 1 2 3 4 6 5  > output_1_0_6_-1_21_-1_2_-2_1_2_3_4_6_5.txt
./matrix_integrate_QCD 1 1 6 -1 21 -1 2 -2 21 1 2 3 4 6 5  > output_1_1_6_-1_21_-1_2_-2_1_2_3_4_6_5.txt
./matrix_integrate_QCD 1 2 6 -1 21 -1 2 -2 21 1 2 3 4 6 5  > output_1_2_6_-1_21_-1_2_-2_1_2_3_4_6_5.txt

./matrix_integrate_QCD 1 0 6 -1 21 -1 2 -2 21 1 3 4 2 6 5  > output_1_0_6_-1_21_-1_2_-2_1_3_4_2_6_5.txt
./matrix_integrate_QCD 1 1 6 -1 21 -1 2 -2 21 1 3 4 2 6 5  > output_1_1_6_-1_21_-1_2_-2_1_3_4_2_6_5.txt
./matrix_integrate_QCD 1 2 6 -1 21 -1 2 -2 21 1 3 4 2 6 5  > output_1_2_6_-1_21_-1_2_-2_1_3_4_2_6_5.txt

# 1

./matrix_integrate_QCD 1 0 6 -1 21 -1 2 -2 21 1 2 6 5 4 3  > output_1_0_6_-1_21_-1_2_-2_1_2_6_5_4_3.txt
./matrix_integrate_QCD 1 1 6 -1 21 -1 2 -2 21 1 2 6 5 4 3  > output_1_1_6_-1_21_-1_2_-2_1_2_6_5_4_3.txt
./matrix_integrate_QCD 1 2 6 -1 21 -1 2 -2 21 1 2 6 5 4 3  > output_1_2_6_-1_21_-1_2_-2_1_2_6_5_4_3.txt

./matrix_integrate_QCD 1 0 6 -1 21 -1 2 -2 21 1 2 5 4 6 3  > output_1_0_6_-1_21_-1_2_-2_1_2_5_4_6_3.txt
./matrix_integrate_QCD 1 1 6 -1 21 -1 2 -2 21 1 2 5 4 6 3  > output_1_1_6_-1_21_-1_2_-2_1_2_5_4_6_3.txt
./matrix_integrate_QCD 1 2 6 -1 21 -1 2 -2 21 1 2 5 4 6 3  > output_1_2_6_-1_21_-1_2_-2_1_2_5_4_6_3.txt

./matrix_integrate_QCD 1 0 6 -1 21 -1 2 -2 21 1 5 4 2 6 3  > output_1_0_6_-1_21_-1_2_-2_1_5_4_2_6_3.txt
./matrix_integrate_QCD 1 1 6 -1 21 -1 2 -2 21 1 5 4 2 6 3  > output_1_1_6_-1_21_-1_2_-2_1_5_4_2_6_3.txt
./matrix_integrate_QCD 1 2 6 -1 21 -1 2 -2 21 1 5 4 2 6 3  > output_1_2_6_-1_21_-1_2_-2_1_5_4_2_6_3.txt


