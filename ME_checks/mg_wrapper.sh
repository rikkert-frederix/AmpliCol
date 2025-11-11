#!/bin/bash

# CHANGE paths of MG folder
path_mg="/home/timea/Documents/MUNKAHELY_LUND/MG5/newest/mg5amcnlo/"

path_code=$(pwd)

version=$(awk '/version/ {print $3}' $path_mg/VERSION)

cp ME_checks/run.sh $path_mg/bin

# Make copy of existing files in MG directory
cp $path_mg/madgraph/iolibs/template_files/matrix_standalone_v4.inc $path_mg/madgraph/iolibs/template_files/matrix_standalone_v4_default.inc
cp $path_mg/madgraph/core/color_algebra.py $path_mg/madgraph/core/color_algebra_default.py

# Copy the modified files depending on version of MG
if [[ "$version" == "3.6.1" ]]; then
    cp ME_checks/matrix_standalone_v4_v361.inc $path_mg/madgraph/iolibs/template_files/matrix_standalone_v4.inc
elif [[ "$version" == "3.6.6" ]]; then
    cp ME_checks/matrix_standalone_v4_v366.inc $path_mg/madgraph/iolibs/template_files/matrix_standalone_v4.inc
else
    echo "Unknown version: $version"
fi

# This is same for the versions
cp color_algebra.py $path_mg/madgraph/core/

# Execute code
cd $path_mg/bin
./run.sh $path_code $@
