#!/usr/bin/sh

# Ensuring we are in the project folder
cd ~/Projets/stage_M2/

# Loading the plots

## The comparison plots
cd comp/
gnuplot comp.plot
cd ../

## The charged system
cd ch-sc/graphs/thermo/
gnuplot thermo.plot

cd ../results/
gnuplot ch-sc.plot

cd ../../../

## The neutral system
cd neutral/graphs/thermo/
gnuplot thermo.plot

cd ../results/
gnuplot neutral.plot

cd ../../../

## The defected system
cd defected/graphs/thermo/
gnuplot thermo.plot

cd ../results/
gnuplot defected.plot

cd ../../../

# Copying the plots to the report folder
cd ~/Projets/stage_M2/

cp -t rapport/src/supercapacitor/ comp/*.pdf
cp -t rapport/src/supercapacitor/ ch-sc/graphs/results/*.pdf ch-sc/graphs/results/1/*.pdf ch-sc/graphs/thermo/*.pdf
cp -t rapport/src/supercapacitor/ neutral/graphs/results/*.pdf neutral/graphs/results/1/*.pdf neutral/graphs/thermo/*.pdf
cp -t rapport/src/supercapacitor/ defected/graphs/results/*.pdf defected/graphs/results/1/*.pdf defected/graphs/thermo/*.pdf
