#!/bin/bash

progress-bar() {
  local duration=${1}


    already_done() { for ((done=0; done<$elapsed; done++)); do printf "########"; done }
    remaining() { for ((remain=$elapsed; remain<$duration; remain++)); do printf "        "; done }
    percentage() { printf "| %s%%" $(( (($elapsed)*100)/($duration)*100/100 )); }
    clean_line() { printf "\r"; }

  for (( elapsed=1; elapsed<=$duration; elapsed++ )); do
      already_done; remaining; percentage
      sleep 1
      clean_line
  done
  clean_line
}

mini-help() {
	echo "
Automated Docking Shell Script using AutoDock Vina
--------------------------------------------------
Created by: June Alexis A. Santos, MSChem Student, AdU

Usage:

First, enable executable permissions by running: chmod +x vinaDock.sh
Then, run the script with by running: bash vinaDock.sh <protein_folder> <ligand_file> <config_file> <alias_enable>

Required flags:

protein_folder - the folder containing one or many protein files which is already prepped for docking.
		 accepted format only in '.pdbqt'
ligand_folder    - the small molecule already prepped
		 accepted format only in '.pdbqt'
config_file    - configuration file containing the BBOX parameters, exhaustiveness, energy diff,
		 and other vina docking options and optimization parameters
alias_enable (-al) - pass vina directory flag if autodock vina was installed using linux based installer

optional flag

--help  -h     - show mini documentation and how-to of this script
"

exit
}

PT=$1
LIG=$2
CONF=$3
AL=$4
WD=$(pwd)

if [[ $AL == "-al" ]]; then
	echo "Fetching prerequisites..."
	echo "This may take awhile..."
	
	sleep 5

	if ! [[ -d $WD/app/ ]]; then
		mkdir $WD/app
		cd $WD/app
		wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.3/vina_1.2.3_linux_x86_64
		
	else
		echo "Prerequisites were already satisfied!"
	fi	

	echo
	echo "Updating app permissions..."
	chmod +x $WD/app/vina_1.2.3_linux_x86_64
	echo "Done!"

	sleep 5
	clear
fi

if [[ $PT == "--help" ]] || [[ $PT == "-h" ]]; then
	mini-help
else
	proteinDir=$WD/$PT
fi

ligand=$WD/$LIG
config=$WD/$CONF

echo "Script: Vina Docking Automation"
echo "Prepared by: June Alexis A. Santos, MS Chemistry Student, AdU"
echo "-------------------------------------------------------------"

echo "Take note: This script will only work with preprocessed ligand
and protein file with specified grid box parameters. I recommend using 
Autodock tools to properly set GridBox params.
-----------------------------------------------------------
 # PDBQT files will only work for ligand and receptor files inputs.
 # This script shall be found on same directory with 
	CONFIG file, protein folder, and Ligand File.

"

sleep 3
echo 
echo "Initializing..."
progress-bar 8
clear

echo
echo "Looking for files..."
sleep 3


#Checking if the files are present. Exiting if otherwise
if [[ -d $proteinDir ]] ; then
	echo 
	echo "FOUND PROTEIN FOLDER: $proteinDir"
else
	echo "Warning: Protein folder not specified. Exiting the script."
	exit
fi
sleep 1

if  [[ -d $ligand ]]; then
	echo
	echo "FOUND LIGAND FILE/S: $ligand"
else
	echo "Warning: Ligand Folder not specified. Exiting the script."
	exit
fi


sleep 3
if  [[ -f $config ]]; then
	echo
	echo "FOUND CONFIG FILE: $config"
	echo 
	echo "Parameters for this run."
	echo "Note: You can always modify this using a text based editor"
	echo "------------------------"
	cat $config
else
	echo "Warning: Config file not specified. Exiting the script."
	exit
fi

sleep 2
echo
echo "------------------------------------"
echo "Creating a target folder container."
echo

if ! [[ -d $WD/Docking_output/ ]]; then
	mkdir $WD/Docking_output
fi

cd $proteinDir


for FILE in find $ligand*; do
	if [[ $FILE == *"pdbqt"* ]]; then
		for x in $find *pdbqt; do
			if [[ $x == *".pdbqt"* ]]; then
				echo "Working on Docking $x"
				sleep 1
				echo "This may take awhile."

				#RUNNING
				if [[ $AL == "-al" ]]; then
					$WD/app/vina_1.2.3_linux_x86_64 --config $config --ligand $FILE --receptor $x --out $WD/Docking_output/${x/".pdbqt"/"_docked.pdbqt"} > $WD/Docking_output/${x/".pdbqt"/""}.txt
				else
					vina --config $config --ligand $FILE --receptor $x --out $WD/Docking_output/${x/".pdbqt"/"_docked.pdbqt"} > $WD/Docking_output/${x/".pdbqt"/""}.txt
				fi

				cat $WD/Docking_output/${x/".pdbqt"/""}.txt

				if ! [[ -f $WD/Docking_output/${x/".pdbqt"/"_docked.pdbqt"} ]]; then
					exit
				else
					echo "Successfully created a pdbqt file output"
				fi
			fi

			echo "Process Done for $x"

		done


	else
			echo "not valid file"
	fi

	echo
done

