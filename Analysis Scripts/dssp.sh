#!/bin/bash

#Run this to get the DSSP data of the protein
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
    if [[ ! -f $directory_name/secondary_structure.xpm || ! -f $directory_name/secondary_structure_2.xpm || ! -f $directory_name/*sec_structure.dat ]]; then
        echo "Running sample: $basename"
        gmx dssp -f "$basename"_mdout_center.xtc -s "$basename"_minimization.tpr -n index.ndx -o "$basename"_sec_structure.dat -xvg none -sel 1

    fi

    cd ..
done

