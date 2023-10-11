import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import glob
import os
from os import path
from tqdm import tqdm
import itertools
import shutil


def get_hbond(traj: str, topol: str, name_code = "", da_distance = 3.2, accepter='whole'):
    # Load trajectory and topology using MDAnalysis
    u = mda.Universe(topol, traj, topology_format='TPR')

    # Define atom selections for donor and acceptor atoms (e.g., based on atom names or residue names)

    acceptor_selection = "resname UNK"   #Ligand resname
    donor_selection = "resid 107:174"   #RBM ResID
    donor_whole = "protein"

    # Calculate hydrogen bonds using MDAnalysis
    # Set cutoff distance to 3.5A (standard) 
    if accepter == 'whole':
        hbonds = HBA(u,
        donors_sel=donor_selection,
        acceptors_sel=acceptor_selection,
        d_a_cutoff=da_distance,
        between=[acceptor_selection,donor_whole],
        d_h_a_angle_cutoff=165,
        update_selections=True)

    elif accepter == 'RBM':
        hbonds = HBA(u, 
        donors_sel=donor_whole, 
        acceptors_sel=acceptor_selection, 
        d_h_a_angle_cutoff=165,
        between=[acceptor_selection,donor_selection],
        d_a_cutoff=da_distance,
        update_selections=True)

    # Initiate HBOND - Analysis
    hbonds.run()

    #Directory and File Mgt
    if not os.path.exists(path.join(path.dirname(topol), f'timestep_hbondcount_{name_code}.csv')):
        # save Hbond with respect to time
        with open(path.join(path.dirname(topol), f'timestep_hbondcount_{name_code}.csv'), 'a') as rd:
            rd.write("time,hbcount\n")
            for ind, count in enumerate(hbonds.count_by_time()):
                rd.write(f"{ind},{count}\n")
            rd.close()

    if not os.path.exists(path.join(path.dirname(topol), f'timtimestep_hbondbyID_{name_code}.txt')):
        # save Hbond unique combination
        with open(path.join(path.dirname(topol), f'timestep_hbondbyID_{name_code}.txt'), 'a') as rd:
            for ind, count_by_id in enumerate(hbonds.count_by_ids()):
                rd.write(f"{ind} - {count_by_id}\n")
            rd.close()

def atomic_groups(traj, topol):
    u = mda.Universe(topol, traj, topology_format='TPR')
    protein = u.atoms.select_atoms("resname UNK")
    protein_with_sol = u.atoms.select_atoms('protein and resid 107:174')
    just_sol = u.atoms.select_atoms("resname CLA and not(around 5.0 protein)")

    for x in list(protein.residues):
        print(x)

def copy_and_rename(directory: str, destination: str):
    if os.path.exists(directory) and os.path.isdir(directory):
        hbonds = glob.glob(os.path.join(directory,"*","*hbondcount*"))
        for file in tqdm(hbonds):
            parent_folder_name = os.path.basename(os.path.dirname(file))
            file_name = os.path.basename(file)
            #copy the file 
            shutil.copy(file,destination)
            #rename the file
            os.rename(os.path.join(destination,file_name),os.path.join(destination,f"{parent_folder_name}_hbcount.csv"))
        print('File Management Done!')

def copy_and_rename_hba(directory: str, destination: str):
    if os.path.exists(directory) and os.path.isdir(directory):
        hbonds = glob.glob(os.path.join(directory,"*","*hbondbyID*"))
        for file in tqdm(hbonds):
            parent_folder_name = os.path.basename(os.path.dirname(file))
            file_name = os.path.basename(file)
            #copy the file 
            shutil.copy(file,destination)
            #rename the file
            os.rename(os.path.join(destination,file_name),os.path.join(destination,f"{parent_folder_name}_hbondbyID.txt"))
        print('File Management Done!')

#########

#Determine HBOND of Ligand to a Residue of Interest (ROI)
def run_hbond_by_lig_roi(traj: str, topol: str, resid=None, name_code = None, savepath=None,da_distance = 3.2):
    #initialize MDA by loading the Topol and Traj
    u = mda.Universe(topol, traj, topology_format='TPR')

    # Define atom selections for donor and acceptor atoms (e.g., based on atom names or residue names)

    #we substract resid with 332 MDA assigns new resid to each Amino acid
    acceptor_selection = "resname UNK"   #Ligand resname
    donor_selection = f"resid {resid-332}"   #Resid of Interest

    # Calculate hydrogen bonds using MDAnalysis
    # Set cutoff distance to 3.5A (standard) 
    hbonds = HBA(u,
        donors_sel=donor_selection,
        acceptors_sel=acceptor_selection,
        d_a_cutoff=da_distance,
        between=[acceptor_selection,donor_selection],
        d_h_a_angle_cutoff=150,
        update_selections=True)

    # Initiate HBOND - Analysis
    print(f'Running post-simulation supplementary HBOND Analysis: LIG-{resid}')
    hbonds.run()


    #Directory and File Mgt
        # save Hbond with respect to time
    with open(path.join(savepath,f'timestep_hbondcount_{name_code}_LIG-{resid}.csv'), 'a') as rd:
        rd.write("time,hbcount\n")
        for ind, count in enumerate(hbonds.count_by_time()):
            rd.write(f"{ind},{count}\n")
        rd.close()

if __name__ == "__main__":
    

    #designed for VM computing. Change directory and conditional designs based on directory and filenaming pointers
    dirpath = '/root/documents/proxim_data'

    #location of post sim hba
    post_hba = "/root/documents/post_hb_res"

    # for contents in os.listdir(dirpath):
    #     if os.path.isdir(path.join(dirpath,contents)):
    #         tpr_file = glob.glob(path.join(dirpath,contents,"*minimization.tpr"))
    #         xtc_file = glob.glob(path.join(dirpath,contents,"*center.xtc"))

    #         get_hbond(xtc_file[0], tpr_file[0], name_code='strong_hb',da_distance=3.5)
    #         # atomic_groups(xtc_file[0], tpr_file[0])
    #         print(f'Finished processing: {path.basename(tpr_file[0])} and {path.basename(xtc_file[0])}')

    #POST SIMULATION RUN
    #set residues of interest
    ROI = {
        "OMIC": [470,492,493,494,452,490,449],
        "ALPHA": [505, 501, 502, 496, 403, 498, 453, 500],
        "N440K": [472, 493, 492, 470, 484, 494, 452, 449, 490],
        "Q493R": [493,456,489,455,494,449,496,453],
        "Q498R": [494,453,497,489,403,498,501,456,455,493,505,495, 496],
        "T478K": [455,502,494,453,449,500,493,497,498,403,501,495,496,505]
        }
    for contents in os.listdir(dirpath):
        if os.path.isdir(path.join(dirpath,contents)):
            tpr_file = glob.glob(path.join(dirpath,contents,"*minimization.tpr"))
            xtc_file = glob.glob(path.join(dirpath,contents,"*center.xtc"))
            if contents in ROI.keys():
                for res in ROI[contents]:
                    run_hbond_by_lig_roi(xtc_file[0], tpr_file[0], resid=res, name_code=contents,savepath=post_hba,da_distance=3.5)
                    print(f'Finished processing: {path.basename(tpr_file[0])} and {path.basename(xtc_file[0])}')
