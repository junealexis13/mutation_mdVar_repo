#!/bin/bash

#custom bash script to calculate per frame distance of Ligand relative to binding site
# Set the folder you want to iterate through
folder=$1
for directory in "$folder"/*/
do
    # Extract the directory name from the path
    directory_name="${directory%/}"

    # Do something with the directory
    echo "Processing directory: $directory_name"
    cd "$directory_name"

    ###CODE###
    current_dir=$(pwd)
    basename=$(basename "$current_dir")

   if [[ ! -f *distance* ]]; then
       echo "Computing DistanceL $basename"
       gmx distance -s "${basename}_minimization.tpr" -f "${basename}_mdout_center.xtc" -select 'resname "UNK"' -select 'resid 436-506' -oav "${basename}_avedistance.xvg"
   fi
    cd ..
done
