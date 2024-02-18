#!/bin/bash

echo "
SLURM SUBMISSION SCRIPT (PT1)
-----------------------------
Author: June Alexis A. Santos, MSc Student, AdU
"

ROOT=$(pwd)

for DIR in */; do
    DIR=$(realpath "$DIR")
    if [[ -d "$DIR" ]]; then
        echo "Entering directory: $DIR"
	cd $DIR
	if ! [[ -f index_unk.ndx ]] ; then
		gmx make_ndx -f $(pwd)/unk.gro -o $(pwd)/index_unk.ndx <<< "0 & ! a H*"$'\n'"q"
        fi

	if ! [[ -f posre_unk.itp ]] ; then
		gmx genrestr -f $(pwd)/unk.gro -n $(pwd)/index_unk.ndx -o $(pwd)/posre_unk.itp -fc 1000 1000 1000 <<< "3"
	fi

        cd ..
    fi
done
