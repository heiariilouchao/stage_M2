#!/usr/bin/sh

# Initializing some variable
project_path='~/Projets/stage_M2/'

# Ensuring we are in the project folder
cd $project_path

# Loading the plots

## The comparison plots
cd comp/
gnuplot comp.plot
cd $project_path

## The charged system
cd ch-sc/graphs/thermo/
gnuplot thermo.plot

cd ../results/
gnuplot ch-sc.plot

## The neutral system
cd neutral/graphs/thermo/
gnuplot thermo.plot

cd ../results/
gnuplot neutral.plot

## The defected system
cd defected/graphs/thermo/
gnuplot thermo.plot

cd ../results/
gnuplot defected.plot

# Copying the plots to the report folder
cd $project_path
cp -t rapport/src/supercapacitor/ comp/*.pdf
cp -t rapport/src/supercapacitor/ ch-sc/graphs/results/*.pdf ch-sc/graphs/results/1/*.pdf ch-sc/graphs/thermo/*.pdf
cp -t rapport/src/supercapacitor/ neutral/graphs/results/*.pdf neutral/graphs/results/1/*.pdf neutral/graphs/thermo/*.pdf
cp -t rapport/src/supercapacitor/ defected/graphs/results/*.pdf defected/graphs/results/1/*.pdf defected/graphs/thermo/*.pdf
