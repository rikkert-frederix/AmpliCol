#! /usr/bin/env sh

make matrix_reweight_QCD

./matrix_reweight_QCD 7 21 21 21 21 21 21 21 1 2 3 4 5 6 7
./matrix_reweight_QCD 7 21 21 21 21 21 21 21 1 3 2 4 5 6 7
./matrix_reweight_QCD 7 21 21 21 21 21 21 21 1 3 4 2 5 6 7
#./matrix_reweight_QCD 7 21 21 21 21 21 21 21 1 3 4 5 2 6 7
#./matrix_reweight_QCD 7 21 21 21 21 21 21 21 1 3 4 5 6 2 7
#./matrix_reweight_QCD 7 21 21 21 21 21 21 21 1 3 4 5 6 7 2




