#!/bin/bash

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
	
	echo "$basename"_minimization.tpr
	echo "Running MMPBSA/decomp assesment on $current_dir"
	gmx_MMPBSA -O -i mmpbsa.in -cs "$basename"_minimization.tpr -ct "$basename"_mdout_center.xtc -ci index.ndx -cg 1 13 -cp topol.top -o mmpbsa_"$basename".dat -eo mmpbsa_"$basename".csv
	cd ..
	
done