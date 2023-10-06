import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.analysis.distances import dist
from MDAnalysis.analysis.contacts import Contacts
import mdtraj as md
import glob, os
from os import path
from tqdm import tqdm
import itertools
import shutil
import numpy as np
import pandas as pd
import time

def get_near_res(traj: str, topol: str, file_code = None):
    # Load trajectory and topology using MDAnalysis
    u = mda.Universe(topol, traj, topology_format='TPR')
    #defining distance cutoff
    distance_cutoff = 2.0 #recommended default for All-Atom MD Simulation by MDA
    
    lig_selection = 'resname UNK'
    rec_selection = 'protein and not (resname UNK)'

    ligand_atoms = u.select_atoms(lig_selection)
    protein_atoms = u.select_atoms(rec_selection)
    #side chain selector
    b_chain = '(not name H*) and (not type NH1) and (not type O) and (not name C)'
    #create a dictionary of frames
    cont = dict()

    a = u.select_atoms(f"(around {distance_cutoff} resname UNK) and (not resname SOL) and (not resname SOD) and (not resname CLA)", updating=True)
    for ts in u.trajectory:
        print(ts.time/1000, list(set([y+332 for y in a.resids])))
    # lig = ligand_atoms.select_atoms("all and not name H*")
    # aa = protein_atoms.select_atoms(f"resid {res.resid}")
    # aa_terminal = aa.select_atoms(f"all and {b_chain}")


    # ndf = pd.DataFrame(cont, index = [x for x in range(1001)])
    # ndf.to_csv(f"others/proximity_check/data_{file_code}.csv")

    # print(get_contacts.results.timeseries[:,1])
    # lres = ligand_atoms.residues
    # rres = protein_atoms.residues

if __name__ == "__main__":
    dirpath = '/Users/junealexissantos/Desktop/PRELIM_DATA'
    save_loc = "/Users/junealexissantos/Desktop/hbnum_Mda/contacts_Mda"

    for contents in os.listdir(dirpath):
        if os.path.isdir(path.join(dirpath,contents)) and "S477N" in contents:
            tpr_file = glob.glob(path.join(dirpath,contents,"*minimization.tpr"))
            xtc_file = glob.glob(path.join(dirpath,contents,"*center.xtc"))
            print(f'Processing {os.path.basename(tpr_file[0])}')
            print(f'Processing {os.path.basename(xtc_file[0])}\n\n')
            time.sleep(2)
            file_name = os.path.basename(path.join(dirpath, contents)).split("_")
            if file_name[-1] == "complex":
                file_code = file_name[0]
            else:
                file_code = file_name[-1]

            get_near_res(xtc_file[0],tpr_file[0], file_code=file_code)


    #Qt = contact per frame / total contacts 
    #FREQUENCY BASED CONTACT ANALYSIS