#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

cd constant/polyMesh
rm -rf boundary points faces neighbour owner
cd ../..
rm -rf data/*.xy
rm -rf data/p.dat
rm -rf data/*.pdf
rm -rf processor*
rm -rf *[1-9]*
rm -rf probes 
rm -rf sets
rm -rf postProcessing
rm log.*
rm snapshot.*
