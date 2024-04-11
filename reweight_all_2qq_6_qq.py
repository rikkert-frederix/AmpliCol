#! /usr/bin/env sh

make matrix_reweight_QCD

./matrix_reweight_QCD 6 -1 1 2 -2 21 21 1 5 6 4 3 2   
./matrix_reweight_QCD 6 -1 1 2 -2 21 21 1 5 4 3 6 2
./matrix_reweight_QCD 6 -1 1 2 -2 21 21 1 4 3 5 6 2

./matrix_reweight_QCD 6 -1 1 2 -2 21 21 1 5 6 2 3 4 
./matrix_reweight_QCD 6 -1 1 2 -2 21 21 1 5 2 3 6 4
./matrix_reweight_QCD 6 -1 1 2 -2 21 21 1 2 3 5 6 4




