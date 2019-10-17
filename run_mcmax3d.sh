#!/bin/bash
cd  /data/users/bportilla/runs/my_MCMax3D_run001
./runmod.sh 
rm *~
N_files=$(ls ../PDS70_runs | wc -l)
cp -r . ../PDS70_runs/run_$(($N_files+1))



