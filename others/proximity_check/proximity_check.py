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

def get_near_res(traj: str, topol: str):
    # Load trajectory and topology using MDAnalysis
    u = mda.Universe(topol, traj, topology_format='TPR')
    #defining distance cutoff
    distance_cutoff = 4.5 #recommended default for All-Atom MD Simulation by MDA
    
    lig_selection = 'resname UNK'
    rec_selection = 'protein and not (resname UNK)'
    ligand_atoms = u.select_atoms(lig_selection)
    protein_atoms = u.select_atoms(rec_selection)

    #create a dictionary of frames
    cont = dict()

    for res in protein_atoms.residues:
        get_contacts = Contacts(u,(lig_selection,f"resname {res.resname}"),refgroup=(ligand_atoms,u.select_atoms(f"resname {res.resname}")),radius=4.5)
        get_contacts.run()
        cont[res.resname+str(res.resid + 332)]= get_contacts.results.timeseries[:,1]

    return pd.DataFrame(cont, [x for x in range(1002)])

    # print(get_contacts.results.timeseries[:,1])
    # lres = ligand_atoms.residues
    # rres = protein_atoms.residues


if __name__ == "__main__":
    dirpath = '/Users/junealexissantos/Desktop/PRELIM_DATA'
    save_loc = "/Users/junealexissantos/Desktop/hbnum_Mda/contacts_Mda"

    for contents in os.listdir(dirpath):
        if os.path.isdir(path.join(dirpath,contents)):
            tpr_file = glob.glob(path.join(dirpath,contents,"*minimization.tpr"))
            xtc_file = glob.glob(path.join(dirpath,contents,"*center.xtc"))


            print(get_near_res(xtc_file[0],tpr_file[0]))

            break