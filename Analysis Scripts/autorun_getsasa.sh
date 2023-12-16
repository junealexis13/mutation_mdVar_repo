#!/bin/bash

#automatically get the SASA from MD trajectories within one folder
# Set the folder you want to iterate through
folder=$1

for directory in "$folder"/*/
do
    # Extract the directory name from the path
    directory_name="${directory%/}"

    # Do something with the directory
    echo "Processing directory: $directory_name"
	cd $directory_name
	
	###CODE###
	current_dir=$(pwd)
	basename=$(basename $current_dir)
	
	if [ -f sasa.xvg ]; then
		echo "Found result file :"$basename"_sasa.csv, skipping MMGBSA runtime"
	else
		echo "$basename"_minimization.tpr
		echo "Running SOLVENT ACCESSIBLE SURFACE AREA assesment on $current_dir"
		echo  "19" | gmx sasa -f "$basename"_mdout_center.xtc -s "$basename"_minimization.tpr -n index.ndx -o "$basename"_sasa.xvg -odg "$basename"_solvation_free_energy.xvg -or "$basename"_avg_perres_sasa.xvg -tv "$basename"_total_v.xvg
	fi
	cd ..
	
done