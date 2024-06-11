#! /usr/bin/env sh

make matrix_reweight_QCD

./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 1 2 4 5 6
./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 2 1 4 5 6
./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 1 4 5 2 6
./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 2 4 5 1 6
./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 4 5 1 2 6
./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 4 5 2 1 6

./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 1 2 6 5 4
./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 2 1 6 5 4
./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 1 6 5 2 4
./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 2 6 5 1 4
./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 6 5 1 2 4
./matrix_reweight_QCD 6 21 21 1 -1 2 -2 3 6 5 2 1 4




