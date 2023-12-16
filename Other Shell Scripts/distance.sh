#!/bin/bash

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
    if [[ ! -f "${basename}_proteinOnlyRMSD.xvg" ]]; then
        echo "Running sample: $basename"
        gmx rms -s "./${basename}_minimization.tpr" -f "./${basename}_mdout_center.xtc" -n ./index.ndx -o "./${basename}_proteinOnlyRMSD.xvg" <<< "4"$'\n'"1"
    fi

    if [[ ! -f "${basename}_calphaRMSD.xvg" ]]; then
        echo "Running sample: $basename"
        gmx rms -s "./${basename}_minimization.tpr" -f "./${basename}_mdout_center.xtc" -n ./index.ndx -o "./${basename}_calphaRMSD.xvg" <<< "4"$'\n'"3"
    fi

#    if [[ ! -f *distance* ]]; then
#        echo "Computing DistanceL $basename"
#        gmx distance -s "${basename}_minimization.tpr" -f "${basename}_mdout_center.xtc" -select 'resname "UNK"' -select 'resid 436-506' -oav "${basename}_avedistance.xvg"
#    fi
    cd ..
done
