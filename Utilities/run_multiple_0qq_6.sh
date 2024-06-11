#! /usr/bin/env sh

make matrix_integrate_QCD

./matrix_integrate_QCD 2 0 6 21 21 21 21 21 21 1 2 3 4 5 6 > output_1_0_6_21_21_21_21_21_21_1_2_3_4_5_6.txt
./matrix_integrate_QCD 2 1 6 21 21 21 21 21 21 1 2 3 4 5 6 > output_1_1_6_21_21_21_21_21_21_1_2_3_4_5_6.txt
./matrix_integrate_QCD 2 2 6 21 21 21 21 21 21 1 2 3 4 5 6 > output_1_2_6_21_21_21_21_21_21_1_2_3_4_5_6.txt

./matrix_integrate_QCD 2 0 6 21 21 21 21 21 21 1 3 2 4 5 6 > output_1_0_6_21_21_21_21_21_21_1_3_2_4_5_6.txt
./matrix_integrate_QCD 2 1 6 21 21 21 21 21 21 1 3 2 4 5 6 > output_1_1_6_21_21_21_21_21_21_1_3_2_4_5_6.txt
./matrix_integrate_QCD 2 2 6 21 21 21 21 21 21 1 3 2 4 5 6 > output_1_2_6_21_21_21_21_21_21_1_3_2_4_5_6.txt

./matrix_integrate_QCD 2 0 6 21 21 21 21 21 21 1 3 4 2 5 6 > output_1_0_6_21_21_21_21_21_21_1_3_4_2_5_6.txt
./matrix_integrate_QCD 2 1 6 21 21 21 21 21 21 1 3 4 2 5 6 > output_1_1_6_21_21_21_21_21_21_1_3_4_2_5_6.txt
./matrix_integrate_QCD 2 2 6 21 21 21 21 21 21 1 3 4 2 5 6 > output_1_2_6_21_21_21_21_21_21_1_3_4_2_5_6.txt

#./matrix_integrate_QCD 2 0 6 21 21 21 21 21 21 1 3 4 5 2 6 > output_1_0_6_21_21_21_21_21_21_1_3_4_5_2_6.txt
#./matrix_integrate_QCD 2 1 6 21 21 21 21 21 21 1 3 4 5 2 6 > output_1_1_6_21_21_21_21_21_21_1_3_4_5_2_6.txt
#./matrix_integrate_QCD 2 2 6 21 21 21 21 21 21 1 3 4 5 2 6 > output_1_2_6_21_21_21_21_21_21_1_3_4_5_2_6.txt

#./matrix_integrate_QCD 2 0 6 21 21 21 21 21 21 1 3 4 5 6 2 > output_1_0_6_21_21_21_21_21_21_1_3_4_5_6_2.txt
#./matrix_integrate_QCD 2 1 6 21 21 21 21 21 21 1 3 4 5 6 2 > output_1_1_6_21_21_21_21_21_21_1_3_4_5_6_2.txt
#./matrix_integrate_QCD 2 2 6 21 21 21 21 21 21 1 3 4 5 6 2 > output_1_2_6_21_21_21_21_21_21_1_3_4_5_6_2.txt
