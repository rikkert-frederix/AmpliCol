#! /usr/bin/env sh

file=$1

cwd='/home/timea/Documents/Uppsala_MUNKAHELY/Projects/Colour_implement/IntegrateGluons/Utilities'
cd "$cwd/compute_unweighting_efficiency" 
./compute_unw_eff $1
