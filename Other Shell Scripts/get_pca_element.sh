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
    if [[ ! -f "$basename"_gibbs.xpm ]]; then
        echo "Running sample: $basename"
        gmx covar -f "$basename"_mdout_center.xtc -s "$basename"_minimization.tpr -n index.ndx -o "$basename"_eigenvl.xvg -tu ns <<< "3"$'\n'"13"
        gmx anaeig -v eigenvec.trr -f "$basename"_mdout_center.xtc -s "$basename"_minimization.tpr -n index.ndx -comp "$basename"_eigcomp.xvg -rmsf "$basename"_eigrmsf.xvg -proj "$basename"_eigproj.xvg -2d "$basename"_2dproj.xvg -tu ns -3d -first 1 -last 3 <<< "3"$'\n'"13"
        gmx sham -f "$basename"_2dproj.xvg -ls "$basename"_gibbs.xpm -lss "$basename"_entropy.xpm -lsh "$basename"_enthalpy.xpm -notime
        ## for Protein
        
        gmx covar -f "$basename"_mdout_center.xtc -s "$basename"_minimization.tpr -n index.ndx -o "$basename"_eigenvl_prot.xvg -v eigenvec_prot.trr -av average_prot.pdb -l covar_prot.log -tu ns <<< "3"$'\n'"1"
        gmx anaeig -v eigenvec_prot.trr -f "$basename"_mdout_center.xtc -s "$basename"_minimization.tpr -n index.ndx -comp "$basename"_eigcomp_prot.xvg -rmsf "$basename"_eigrmsf_prot.xvg -proj "$basename"_eigproj_prot.xvg -2d "$basename"_2dproj_prot.xvg -tu ns -3d 3dproj_prot.pdb -first 1 -last 3 <<< "3"$'\n'"1"
        gmx sham -f "$basename"_2dproj_prot.xvg -ls "$basename"_gibbs_prot.xpm -lss "$basename"_entropy_prot.xpm -lsh "$basename"_enthalpy_prot.xpm -notime
    else
        echo "Skipping this file. FEL found."
    fi

    cd ..
    
echo "removing duplicate backups"
rm */#*#
echo "PROCESS DONE! $(date)"

done


