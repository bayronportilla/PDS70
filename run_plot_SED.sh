#!/bin/bash

python3 plot_SED.py $1 $2 $3
mkdir plots/run_$1_$2_$3
mv fig.png plots/run_$1_$2_$3/run_$1_$2_$3.png
for item in "$@"; 
do    
    rm spectrum_PDS70_star_$item.dat
    rm spectrum_PDS70_system_$item.dat
done
rm *~

