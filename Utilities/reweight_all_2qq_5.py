#! /usr/bin/env sh

make matrix_reweight_QCD

./matrix_reweight_QCD 5 -1 1 2 -2 21 1 5 4 3 2   
./matrix_reweight_QCD 5 -1 1 2 -2 21 1 4 3 5 2

./matrix_reweight_QCD 5 -1 1 2 -2 21 1 5 2 3 4
./matrix_reweight_QCD 5 -1 1 2 -2 21 1 2 3 5 4



