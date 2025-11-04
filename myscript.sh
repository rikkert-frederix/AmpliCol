#!/bin/bash

# CHANGE paths of MG folder and code 
path_mg="/home/timea/Documents/MUNKAHELY_LUND/MG5/LC_MG_v3_6_1/"
path_code="/home/timea/Documents/Uppsala_MUNKAHELY/Projects/Colour_implement/new_integrator/IntegrateGluons/"

cp run.sh $path_mg/bin
cp matrix_standalone_v4.inc $path_mg/madgraph/iolibs/template_files
cp color_algebra.py $path_mg/madgraph/core/
cd $path_mg/bin

./run.sh $path_code $@
