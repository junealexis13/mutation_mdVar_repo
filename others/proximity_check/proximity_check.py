import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.analysis.distances import dist
from MDAnalysis.analysis.contacts import Contacts
import matplotlib.pyplot as plt
import glob, os
from os import path
from tqdm import tqdm
import itertools
import shutil
import numpy as np
import pandas as pd
import time

def get_near_res(traj: str, topol: str, savepath=None, file_code = None):
    # Load trajectory and topology using MDAnalysis
    u = mda.Universe(topol, traj, topology_format='TPR')
    #defining distance cutoff
    distance_cutoff = 3.0 #recommended default for All-Atom MD Simulation by MDA
    
    lig_selection = 'resname UNK'
    rec_selection = 'protein and not (resname UNK)'

    ligand_atoms = u.select_atoms(lig_selection)
    protein_atoms = u.select_atoms(rec_selection)
    #side chain selector
    b_chain = '(not name H*) and (not type NH1) and (not type O) and (not name C)'

    #Creating a trajectory cutoff | This will halt the proximity calculation at specified frame
    mut = {
        'N501Y': 170,
        'E484A': 102,
        'Q493R': 626,
        'Q498R': 611,
        'T478K': 373,
        'N440K': 508,
        'Y505H': 385,
        'ALPH' : 75
    }
    #create a dictionary of frames
    cont = dict()

    
    start_timer = time.perf_counter()

    a = u.select_atoms(f"(around {distance_cutoff} resname UNK) and (not resname SOL) and (not resname SOD) and (not resname CLA)", updating=True)

    #create a file to save each frame data
    with open(f"{savepath}/{file_code}.csv","a") as rd:

        if file_code.upper() not in mut.keys():
            for ts in tqdm(u.trajectory, desc='Loading Traj.'):
                cont[ts.time/100] = list(set([f"{y+332}{x}" for y,x in zip(a.resids,a.resnames)]))
                rd.write(",".join([str(ts.time/100)]+ list(set([f"{y+332}{x}" for y,x in zip(a.resids,a.resnames)]))))
                rd.write("\n")
        elif file_code.upper() in mut.keys():
            print(f"Note: {file_code.upper()} has an assigned breakpoint...")   
            for ts in tqdm(u.trajectory[:mut[file_code.upper()]]):
                cont[ts.time/100] = list(set([f"{y+332}{x}" for y,x in zip(a.resids,a.resnames)]))
                rd.write(",".join([str(ts.time/100)]+ list(set([f"{y+332}{x}" for y,x in zip(a.resids,a.resnames)]))))
                rd.write("\n")
        rd.close()


    end_timer = time.perf_counter()
    
    #Printouts
    print(60*"-")
    print(f"Finished Loading: took {end_timer - start_timer:.2f}seconds to load.")
    return cont


#################

def process_data(data: dict):
    print('Fetching contact distribution...')
    time.sleep(0.5)
    #store unique interacting AA
    interacting_AA = []

    #store each frame
    frame_AA = {}

    #Set timer again
    start_timer = time.perf_counter()

    for frame, res in tqdm(data.items()):
        #unique AA search
        [interacting_AA.append(x) for x in sorted(res) if x not in interacting_AA]

    print("Processing data... This may take a while.")

    for frame, res in data.items():
        #matching with unique AA and this could be slow 
        for ires in sorted(interacting_AA):
            if ires in res:
                if ires not in frame_AA.keys():
                    frame_AA[ires] = 1
                else:
                    frame_AA[ires] += 1
    end_timer = time.perf_counter()
    print(f"Finished Loading: took {end_timer - start_timer:.2f}seconds to load.")
    return frame_AA

#################

def plot_aa(data_ff: dict, sample_code: str, savepath=None, show=False):

    #transform into df
    df = pd.DataFrame.from_dict(data_ff, orient='index')
    df.columns = ['count']
    df_sorted = df.sort_index()

    fig, ax = plt.subplots(figsize=(16,8))
    ax.bar(df_sorted.index, df_sorted['count'])
    ax.tick_params(which="major", axis="x", labelrotation=90)
    ax.set_title(f"Distribution of interacting amino acid with ligand for {sample_code} complex")
    ax.set_ylim(0, 1000)
    plt.savefig(f"{savepath}/{file_code}.jpg",dpi=180)

    if show:
        print("Viewing plot...")
        plt.show()




if __name__ == "__main__":
    dirpath = '/Users/junealexissantos/Desktop/PRELIM_DATA'
    save_loc = "/Users/junealexissantos/Desktop/savepath_interacting_aa/"

    for contents in os.listdir(dirpath):
        if os.path.isdir(path.join(dirpath,contents)):
            tpr_file = glob.glob(path.join(dirpath,contents,"*minimization.tpr"))
            xtc_file = glob.glob(path.join(dirpath,contents,"*center.xtc"))
            print(120*"#")
            print(f'\nProcessing {os.path.basename(tpr_file[0])}')
            print(f'Processing {os.path.basename(xtc_file[0])}\n')
            time.sleep(2)
            file_name = os.path.basename(path.join(dirpath, contents)).split("_")
            if file_name[-1] == "complex":
                file_code = file_name[0]
            else:
                file_code = file_name[-1]

            data = get_near_res(xtc_file[0],tpr_file[0], savepath=save_loc, file_code=file_code)
            aa_distribution = process_data(data)
            plot_aa(aa_distribution, sample_code=file_code, savepath=save_loc)


            


    #Qt = contact per frame / total contacts 
    #FREQUENCY BASED CONTACT ANALYSIS