#!/bin/sh


OUTPUT="output"
DATA="data"

if [ ! -d $DATA ]
then
    exit 1;
fi

if [ ! -d $OUTPUT ]
then
    mkdir $OUTPUT
fi

source ~/Python/data/bin/activate
for file in $DATA/log.*
do
    base=${file#$DATA/log.}
    echo "Processing $base..."
    python ~/work/stage_M2/utils/extract_thermo.py $file $OUTPUT/$base.log
done