#! /usr/bin/env sh

make 
# 2

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 1 2 4 5 6  > output_1_0_6_21_21_1_-1_2_-2_3_1_2_4_5_6.txt
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 1 2 4 5 6  > output_1_1_6_21_21_1_-1_2_-2_3_1_2_4_5_6.txt
./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 1 2 4 5 6  > output_1_2_6_21_21_1_-1_2_-2_3_1_2_4_5_6.txt

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 2 1 4 5 6  > output_1_0_6_21_21_1_-1_2_-2_3_2_1_4_5_6.txt
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 2 1 4 5 6  > output_1_1_6_21_21_1_-1_2_-2_3_2_1_4_5_6.txt
./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 2 1 4 5 6  > output_1_2_6_21_21_1_-1_2_-2_3_2_1_4_5_6.txt

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 1 4 5 2 6  > output_1_0_6_21_21_1_-1_2_-2_3_1_4_5_2_6.txt
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 1 4 5 2 6  > output_1_1_6_21_21_1_-1_2_-2_3_1_4_5_2_6.txt
./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 1 4 5 2 6  > output_1_2_6_21_21_1_-1_2_-2_3_1_4_5_2_6.txt

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 2 4 5 1 6  > output_1_0_6_21_21_1_-1_2_-2_3_2_4_5_1_6.txt
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 2 4 5 1 6  > output_1_1_6_21_21_1_-1_2_-2_3_2_4_5_1_6.txt
./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 2 4 5 1 6  > output_1_2_6_21_21_1_-1_2_-2_3_2_4_5_1_6.txt

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 4 5 1 2 6  > output_1_0_6_21_21_1_-1_2_-2_3_4_5_1_2_6.txt
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 4 5 1 2 6  > output_1_1_6_21_21_1_-1_2_-2_3_4_5_1_2_6.txt#./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 4 5 1 2 6  > output_1_2_6_21_21_1_-1_2_-2_3_4_5_1_2_6.txt

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 4 5 2 1 6  > output_1_0_6_21_21_1_-1_2_-2_3_4_5_2_1_6.tx
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 4 5 2 1 6  > output_1_1_6_21_21_1_-1_2_-2_3_4_5_2_1_6.txt
./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 4 5 2 1 6  > output_1_2_6_21_21_1_-1_2_-2_3_4_5_2_1_6.txt

# 1

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 1 2 6 5 4  > output_1_0_6_21_21_1_-1_2_-2_3_1_2_6_5_4.txt
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 1 2 6 5 4  > output_1_1_6_21_21_1_-1_2_-2_3_1_2_6_5_4.txt
./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 1 2 6 5 4  > output_1_2_6_21_21_1_-1_2_-2_3_1_2_6_5_4.txt

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 2 1 6 5 4  > output_1_0_6_21_21_1_-1_2_-2_3_1_2_6_5_4.txt
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 2 1 6 5 4  > output_1_1_6_21_21_1_-1_2_-2_3_1_2_6_5_4.txt
./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 2 1 6 5 4  > output_1_2_6_21_21_1_-1_2_-2_3_1_2_6_5_4.txt

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 1 6 5 2 4  > output_1_0_6_21_21_1_-1_2_-2_3_1_6_5_2_4.txt
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 1 6 5 2 4  > output_1_1_6_21_21_1_-1_2_-2_3_1_6_5_2_4.txt
./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 1 6 5 2 4  > output_1_2_6_21_21_1_-1_2_-2_3_1_6_5_124.txt

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 2 6 5 1 4  > output_1_0_6_21_21_1_-1_2_-2_3_2_6_5_1_4.txt
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 2 6 5 1 4  > output_1_1_6_21_21_1_-1_2_-2_3_2_6_5_1_4.txt
./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 2 6 5 1 4  > output_1_2_6_21_21_1_-1_2_-2_3_2_6_5_1_4.txt

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 6 5 1 2 4  > output_1_0_6_21_21_1_-1_2_-2_3_6_5_1_2_4.txt
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 6 5 1 2 4  > output_1_1_6_21_21_1_-1_2_-2_3_6_5_1_2_4.txt
./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 6 5 1 2 4  > output_1_2_6_21_21_1_-1_2_-2_3_6_5_1_2_4.txt

./matrix_integrate_QCD 1 0 6 21 21 1 -1 2 -2 3 6 5 2 1 4  > output_1_0_6_21_21_1_-1_2_-2_3_6_5_2_1_4.txt
./matrix_integrate_QCD 1 1 6 21 21 1 -1 2 -2 3 6 5 2 1 4  > output_1_1_6_21_21_1_-1_2_-2_3_6_5_2_1_4.txt
./matrix_integrate_QCD 1 2 6 21 21 1 -1 2 -2 3 6 5 2 1 4  > output_1_2_6_21_21_1_-1_2_-2_3_6_5_2_1_4.txt


