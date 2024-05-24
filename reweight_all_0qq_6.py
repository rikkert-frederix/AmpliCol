#! /usr/bin/env sh

make matrix_reweight_QCD

./matrix_reweight_QCD 6 21 21 21 21 21 21 1 2 3 4 5 6
./matrix_reweight_QCD 6 21 21 21 21 21 21 1 3 2 4 5 6
./matrix_reweight_QCD 6 21 21 21 21 21 21 1 3 4 2 5 6
./matrix_reweight_QCD 6 21 21 21 21 21 21 1 3 4 5 2 6
./matrix_reweight_QCD 6 21 21 21 21 21 21 1 3 4 5 6 2



