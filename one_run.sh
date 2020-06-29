th#!/bin/bash

cd /data/users/bportilla/runs/final_runs/$1
module load anaconda3
rm sim_log.txt
rm diff.txt
touch sim_log.txt
touch diff.txt
echo $2>>diff.txt
echo "running Monte Carlo">>sim_log.txt
./run_thermal.sh 
echo "Raytracing at ALMA wavelenght">>sim_log.txt
./run_image_alma.sh
echo "Raytracing at SPHERE wavelenght">>sim_log.txt
./run_image_jband.sh 
echo "Running PDS70 pipeline">>sim_log.txt
python3 PDS70_pipeline.py
: '
echo "Running cflux for ALMA">>sim_log.txt
python3 cflux_alma_local.py
echo "Running cflux for SPHERE">>sim_log.txt
python3 cflux_jband_local.py 
echo "Creating image with final results">>sim_log.txt
python3 profiles_modeled_local.py 
echo "Finding ratios for cflux radial profile">>sim_log.txt
python3 jband_ratios_local.py
echo "Preparing alma image">>sim_log.txt
python3 prepare_images_local.py
echo "Extracting azimuthal alma profile">>sim_log.txt
python3 azprofile_local.py
echo "Getting ratios">>sim_log.txt
python3 ratio.py
'
echo "Everything done">>sim_log.txt





