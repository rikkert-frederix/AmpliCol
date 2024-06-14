#! /usr/bin/env sh

file=$1

cd "./compute_unweighting_efficiency" 
pwd
echo $1
./compute_unw_eff $1
