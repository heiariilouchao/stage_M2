!/usr/bin/sh

# Initializing some options
script_opt='--start=4000000'

# Ensuring we are in the project folder
cd ~/Projets/stage_M2/

# Going in the data processing folder
cd traitement_donnees/

## Removing every previous outputs
rm output/graphite/* output/graphite-rdf/* output/graphite-defect/*

## Recompiling the binaries
make reset
make graphite
make graphite-rdf
make defect

# Processing the data

## The charged system
bin/graphite $script_opt data/ch-sc/1.lammpstrj
bin/graphite-rdf $script_opt data/ch-sc/1.lammpstrj

cd ..
mv -t ch-sc/output/1/ traitement_donnees/output/graphite/* traitement_donnees/output/graphite-rdf/*

cd traitement_donnees/
bin/graphite $script_opt data/ch-sc/2.lammpstrj
bin/graphite-rdf $script_opt data/ch-sc/2.lammpstrj

cd ..
mv -t ch-sc/output/2/ traitement_donnees/output/graphite/* traitement_donnees/output/graphite-rdf/*
cd traitement_donnees/

## The neutral system
bin/graphite $script_opt data/neutral/1.lammpstrj
bin/graphite-rdf $script_opt data/neutral/1.lammpstrj

cd ..
mv -t neutral/output/1/ traitement_donnees/output/graphite/* traitement_donnees/output/graphite-rdf/*

cd traitement_donnees/
bin/graphite $script_opt data/neutral/2.lammpstrj
bin/graphite-rdf $script_opt data/neutral/2.lammpstrj

cd ..
mv -t neutral/output/2/ traitement_donnees/output/graphite/* traitement_donnees/output/graphite-rdf/*
cd traitement_donnees/

## The defected system
bin/graphite $script_opt data/defected/1.lammpstrj
bin/graphite-rdf $script_opt data/defected/1.lammpstrj
bin/defect $script_opt data/defected/1.lammpstrj

cd ..
mv -t defected/output/1/ traitement_donnees/output/graphite/* traitement_donnees/output/graphite-rdf/* traitement_donnees/output/graphite-defect/*

cd traitement_donnees/
bin/graphite $script_opt data/defected/2.lammpstrj
bin/graphite-rdf $script_opt data/defected/2.lammpstrj
bin/defect $script_opt data/defected/2.lammpstrj

cd ..
mv -t defected/output/2/ traitement_donnees/output/graphite/* traitement_donnees/output/graphite-rdf/* traitement_donnees/output/graphite-defect/*
