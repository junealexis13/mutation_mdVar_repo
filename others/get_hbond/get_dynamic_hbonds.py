#dynamic hydrogen bonding script - used to monitor the hydrogen bond form and identify the involved atoms and residues
import os, glob
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import numpy as np

class Dynamic_HBOND_Assessment:
    def __init__(self, traj: str, topol: str):
        self.TRAJ = traj
        self.TOPOL = topol
    
    def get_hbond(self, da_distance = 3.2, accepter='whole'):
        # Load trajectory and topology using MDAnalysis
        u = mda.Universe(self.TOPOL, self.TRAJ, topology_format='TPR')

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
        return self.get_atoms(hbonds.count_by_ids())


    def get_atoms(self, count_by_id: np.ndarray ):
        # Load trajectory and topology using MDAnalysis
        u = mda.Universe(self.TOPOL, self.TRAJ, topology_format='TPR')
        dynamic_int_hbonds = {}
        ncounter = 0
        for hbond_mdata in count_by_id:
            bond_interaction = []
            for i, info in enumerate(hbond_mdata):
                # the atom involved
                atom_ = u.atoms[info - 1]
                if i == 0 or i == 1 or i == 2:
                    bond_interaction.append(f"{atom_.resname}{atom_.resid}({atom_.element})")
                elif i == 3:
                    dynamic_int_hbonds[ncounter] = {
                        "bond_data": bond_interaction,
                        "count": info
                    }
            ncounter += 1
        return dynamic_int_hbonds


