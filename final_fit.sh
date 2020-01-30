#!/bin/bash
: '
# Create run dir
rm -rf /data/users/bportilla/runs/run_amax
mkdir /data/users/bportilla/runs/run_amax
mkdir /data/users/bportilla/runs/run_amax/run_amax_initial

# Change to new dir
cd /data/users/bportilla/runs/run_amax/run_amax_initial

# Create log file
touch sim_log.txt
echo "Hello!">>sim_log.txt
echo "log file created">>sim_log.txt

# Load python3 and create density profile
module load anaconda3
echo "Anaconda3 loaded">>sim_log.txt
python3 /Users/users/bportilla/Documents/first_project/scripts/PDS70/surface_density/dprofile.py
echo "Density profile created">>sim_log.txt

# Copy run template (MCMax3D_images)
cp -r /data/users/bportilla/runs/MCMax3D_images/* .    

# Copy aditional staff from PDS 70 dir
cp /Users/users/bportilla/Documents/first_project/scripts/PDS70/alma_radial_profile_observed.dat .
cp /Users/users/bportilla/Documents/first_project/scripts/PDS70/jband_radial_profile_observed.dat .
cp /Users/users/bportilla/Documents/first_project/scripts/PDS70/scripts_paper_01/cflux_alma_local.py .
cp /Users/users/bportilla/Documents/first_project/scripts/PDS70/scripts_paper_01/cflux_jband_local.py .
cp /Users/users/bportilla/Documents/first_project/scripts/PDS70/scripts_paper_01/interpol_local.py .
cp /Users/users/bportilla/Documents/first_project/scripts/PDS70/scripts_paper_01/profiles_modeled_local.py .
cp /Users/users/bportilla/Documents/first_project/scripts/PDS70/PDS70_photometry.dat .
cp /Users/users/bportilla/Documents/first_project/scripts/PDS70/scripts_paper_01/cdm.py .

echo "Doing mother run">>sim_log.txt
# Mother run
./run_thermal.sh 
echo "Monte Carlo run at three level done">>sim_log.txt
./run_image_alma.sh
echo "Alma image raytraced at three level">>sim_log.txt
./run_image_jband.sh
echo "jband image raytraced at three level">>sim_log.txt
    
# Analysis at three level
python3 cflux_alma_local.py
echo "cflux of alma done at three level">>sim_log.txt
python3 cflux_jband_local.py
echo "cflux of jband done at three level">>sim_log.txt

lim_alma_error=28
error_alma=$(python3 cdm.py)
echo "Initial alma error computed">>sim_log.txt
echo "$error_alma">>sim_log.txt
echo "Starting alma fitting at three level">>sim_log.txt

# Fitting alma radial profile at three level
while [ $error_alma -gt $lim_alma_error ]
do
    python3 interpol_local.py 
    ./run_thermal.sh 
    ./run_image_alma.sh
    ./run_image_jband.sh
    python3 cflux_alma_local.py
    python3 cflux_jband_local.py
    current_error=$(python3 cdm.py)
    let error_alma=$current_error
    echo "Iteration finished">>sim_log.txt
    echo "$error_alma">>sim_log.txt
done

echo "Iterations at three level done">>sim_log.txt

# Plotting result at three level
python3 profiles_modeled_local.py
echo "Everything done at three level">>sim_log.txt
'
cd /data/users/bportilla/runs/run_amax

# Creating folders for varying amax 
for var in 10.0 1.0 0.1 0.05
do 
    rm -rf /data/users/bportilla/runs/run_amax/amax_$var
    mkdir /data/users/bportilla/runs/run_amax/amax_$var
    cd /data/users/bportilla/runs/run_amax/amax_$var
    module load anaconda3
    cp -r /data/users/bportilla/runs/run_amax/run_amax_initial/* .
    rm sim_log.txt input.dat 
    touch sim_log_$var.txt
    echo "sim_log_$var.txt created">>sim_log_$var.txt
    cp /Users/users/bportilla/Documents/first_project/scripts/PDS70/scripts_paper_01/gen_input_z1.py .
    python3 gen_input_z1.py $var
    echo "New input file generated for amax=$var">>sim_log_$var.txt
    echo "Running Monte Carlo"
    ./run_thermal.sh
    echo "Monte Carlo done">>sim_log_$var.txt
    ./run_image_alma.sh
    echo "Alma image raytraced">>sim_log_$var.txt
    ./run_image_jband.sh
    echo "jband image raytraced">>sim_log_$var.txt
    # Analysis for new amax
    python3 cflux_alma_local.py
    echo "cflux of alma done">>sim_log_$var.txt
    python3 cflux_jband_local.py
    echo "cflux of jband done">>sim_log_$var.txt
    error_alma=$(python3 cdm.py)
    echo "Initial alma error computed">>sim_log_$var.txt
    echo "$error_alma">>sim_log_$var.txt
    echo "Starting alma fitting at amax=$var">>sim_log_$var.txt
    # Fitting alma radial profile for particular amax
    niter=0
    while [ $error_alma -gt $lim_alma_error ] && [ $niter -lt 4]
    do
	python3 interpol_local.py 
	./run_thermal.sh 
	./run_image_alma.sh
	./run_image_jband.sh
	python3 cflux_alma_local.py
	python3 cflux_jband_local.py
	current_error=$(python3 cdm.py)
	let error_alma=$current_error
	echo "Iteration finished">>sim_log_$var.txt
	echo "$error_alma">>sim_log_$var.txt
	((niter++))
    done
    echo "Iterations done">>sim_log_$var.txt
    python3 profiles_modeled_local.py
    echo "Everything done at amax=$var">>sim_log_$var.txt
done
