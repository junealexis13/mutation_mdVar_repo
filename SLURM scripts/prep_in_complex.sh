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
Script for Preparing Protein-Ligand Complex
--------------------------------------------------
Created by: June Alexis A. Santos, MSChem Student, AdU

Usage:

First, enable executable permissions by running: chmod +x prep_in_complex.sh
Then, run the script with bash by running: bash prep_in_complex.sh <complex_folder> <forcefield_get> <ligand_code> <download_mode>

Required flags:

complex_folder - folder containing all the protein-ligand complex file in pdb format

forcefield_get - option enabling the user to choose between getting forcefield manually or by using CHARMM36 by default
                  (NA) - No Link will be passed. The user will be required to pass a prompt on default gromacs forcefield package.
                  (DEF) - By default, the updated CHARMM36 forcefield will be used on creating a topology compliant file (.gro)
                  (LINK) - Pass a forcefield link and the script will use the file as main Forcefield to be applied.

ligand_code - ligand code to be used as basis for extraction. Should be all caps.

optional flag:

download_mode - option to automatically download required prerequisites.
                (1) Update the system and download pre-requisite resources


--help  -h     - show mini documentation and how-to of this script
"

exit
}

#__init__

DIR=$1
LINK=$2
LIG=$3
DL=$4

ROOT=$(pwd)

if [[ $DIR == "-h" ]] || [[ $DIR == "--help" ]]; then
  mini-help
fi



echo "Script for Preparing Protein-Ligand Complex
#-----------------------------------------------------------#
Prepared by: June Alexis A. Santos, MS Chemistry Student, AdU

Take note: This script will only work with preprocessed Protein-Ligand
complex. The complex file should be presented in .pdb format.
-----------------------------------------------------------

"
sleep 1
echo "Initializing Script"
progress-bar 3

clear

#actual directory containing files to be prepped

echo "The working directory containing all the complexes was set to $ROOT."

sleep 2


#installing dependencies
echo "Attempting to install dependencies..."
echo "--------------------------------------------"
echo "please enter your terminal password."

if [[ $DL == "1" ]]; then
  sudo apt-get update && clear && sudo apt-get install -y gromacs && sudo apt-get install -y openbabel && sudo apt-get install -y perl
  sudo apt install python3
  pip install networkx==2.3
  pip install numpy
fi

clear

echo "PHASE I -- DOWNLOADING NECESSARY FILES FOR MD"
echo "---------------------------------------------"

sleep 2
echo "Step 1.1 - Acquiring Force Field"
echo 
sleep 1

if [ "$LINK" = 'NA' ]; then
	echo "No link was passed"
elif [ "$LINK" = 'DEF' ]; then
	echo "Default Force Field Selected. Downloading LATEST - CHARMM36 forcefield"
	sleep 2
	cd $DIR
	if ! [[ -f default_forcefield.ff.tgz ]] || [[ -d default_forcefield.ff ]]; then
		wget -O default_forcefield.ff.tgz http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz
 	else
		echo "Forcefield file already exist. Skipping the download scheme."
	fi 
	echo "DONE"


else
	echo "Downloading the desired forcefield"
 	echo "Link passed. $LINK"
	sleep 1
	cd $DIR && wget $LINK
	sleep 1
fi

echo "Step 1.2 -- Getting Necesary Scripts"
sleep 2
echo "Getting PERL Sort and CHARMM2GMX Scripts"
echo "Checking availability..."
sleep 2
if ! [[ -f sort_bonds.pl ]]; then
	wget -O sort_bonds.pl http://www.mdtutorials.com/gmx/complex/Files/sort_mol2_bonds.pl
	wget -O cgenff_charmm2gmx_py3_nx2.py https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/cgenff_charmm2gmx_py3_nx2.py
else
	echo "Perl script is present.Skipping..."
fi

echo "Step 2 - Directory Management"
echo "-----------------------------"
sleep 1
echo "Sorting the files and creating a custom folder for each samples"
sleep 1
echo "Note: This script assumes that the ligand complex was in PDB format.
If files were not yet in proper format, CANCEL This Script by pressing CTRL + C or CMD + C."

echo "-----------------------------"
echo "Sorting the directory"
sleep 1

echo "Finding PDB files..."
sleep 1.5

for files in $(find *.pdb); do
	fName=${files/".pdb"/""}
	echo "FOUND: $files" && mkdir $fName
  	echo "-----------------------------"
  	echo "Created new WSPACE: $DIR$fName/"
  	sleep 1
  	echo "Moving the Force Field FILE and PDBFILE to each new working Dir: -->  $DIR/$fName..."
	

	if [[ -d default_forcefield.ff ]]; then
		cp -r default_forcefield.ff/ $DIR/$fName/
		cp cgenff_charmm2gmx_py3_nx2.py $DIR/$fName  	
	elif [[ -f $(find *ff.tgz) ]]; then
		mv $files $DIR/$fName && cp $(find *ff.tgz) $DIR/$fName
		pushd $fName
  		echo "Unpacking ForceField..."
  		tar -xvzf $(find *ff.tgz)
		cp ../cgenff_charmm2gmx_py3_nx2.py ../$fName
  		rm $(find *ff.tgz)
	else
		echo "Error Fetching FF files"
		exit
	fi

  	echo "Step 3 - Creating a Force Field Compliant Coord. File"
  	echo "-----------------------------"
  	sleep 1 && clear

  	sleep 1
  	echo "Attempting to isolate the NON-FF Ligand from the complex."
  	echo "-----------------------------"
  	sleep 1
  	echo "Processing..."
  	echo
  	echo "Removing non bonding HOH"
  	grep -v HOH $files > "${files/".pdb"/"_cleaned.pdb"}" && grep -v $LIG $(find *_cleaned.pdb) > "${files/".pdb"/"_prot.pdb"}"
  	sleep 2
  	echo "Isolating the Ligand..."
  	grep $LIG $(find *_cleaned.pdb) > ${files/".pdb"/"_"}$LIG.pdb    #Ligand File
  	echo "-----------------------------"
  	echo "Writing Protein Topology using the available Force Field."
  	echo "-----------------------------"
  	if [[ $LINK == "DEF" ]]; then
    		pfile=$(find *_prot.pdb)
    		gmx pdb2gmx -f $pfile -ff charmm36-jul2022 -water tip3p -o $(pwd)/$fName.gro -ignh
  	else
    		pfile=$(find *_prot.pdb)
    		ffFile=$(find . -name *.ff | cut -c 3-)
    		ffTag=$(ffFile/".ff"/"")
    		gmx pdb2gmx -f $pfile -ff $ffTag -water tip3p -o $(pwd)/$fName.gro -ignh
  	fi


  	echo "Post MOL2 Processing"
  	sleep 2
	LIGNAME=${files/".pdb"/"_"}$LIG
  	obabel $LIGNAME.pdb -O $LIGNAME.mol2 -h			#this converts and addH
  	sed "s/$LIGNAME.pdb/$LIG/g" $LIGNAME.mol2 > $LIG.mol2
	echo "DONE!"
	mv $files ../
  	
	cp ../sort_bonds.pl .
	sleep 2
	if [[ -f $LIG.mol2 ]]; then
		echo "Sorting bonds..."
		perl sort_bonds.pl $LIG.mol2 sorted.mol2
	fi

	popd 

done

echo "Processing done"

