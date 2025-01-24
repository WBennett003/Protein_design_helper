import os, sys

from pyrosetta import *
from pyrosetta.rosetta import *

parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(parent, 'include'))
from silent_tools import silent_tools

init( "-beta_nov16 -in:file:silent_struct_type binary -mute all" +
    " -use_terminal_residues true -mute basic.io.database core.scoring" )

import math
import numpy as np
from collections import defaultdict, OrderedDict
from scipy.spatial import cKDTree
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, DSSP, Selection, Polypeptide, PDBIO, Select, Chain, Superimposer
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.Polypeptide import is_aa

import pathlib
import argparse
import json

three_to_one_map = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# identify interacting residues at the binder interface
def hotspot_residues(trajectory_pdb, binder_chain="B", atom_distance_cutoff=4.0):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", trajectory_pdb)

    # Get the specified chain
    binder_atoms = Selection.unfold_entities(structure[0][binder_chain], 'A')
    binder_coords = np.array([atom.coord for atom in binder_atoms])

    # Get atoms and coords for the target chain
    target_atoms = Selection.unfold_entities(structure[0]['A'], 'A')
    target_coords = np.array([atom.coord for atom in target_atoms])

    # Build KD trees for both chains
    binder_tree = cKDTree(binder_coords)
    target_tree = cKDTree(target_coords)

    # Prepare to collect interacting residues
    interacting_residues = {}

    # Query the tree for pairs of atoms within the distance cutoff
    pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)

    # Process each binder atom's interactions
    for binder_idx, close_indices in enumerate(pairs):
        binder_residue = binder_atoms[binder_idx].get_parent()
        binder_resname = binder_residue.get_resname()

        # Convert three-letter code to single-letter code using the manual dictionary
        if binder_resname in three_to_one_map:
            aa_single_letter = three_to_one_map[binder_resname]
            for close_idx in close_indices:
                target_residue = target_atoms[close_idx].get_parent()
                interacting_residues[binder_residue.id[1]] = aa_single_letter

    return interacting_residues


class StructManager():
    '''
    This class handles all of the input and output for the ProteinMPNN model. It deals with silent files vs. pdbs,
    checkpointing, and writing of outputs

    Note: This class could be moved to a separate file
    '''

    def __init__(self, args):
        self.args = args

        self.silent = False
        if not args.silent == '':
            self.silent = True

            self.struct_iterator = silent_tools.get_silent_index(args.silent)['tags']

            self.sfd_in = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
            self.sfd_in.read_file(args.silent)

            self.sfd_out = core.io.silent.SilentFileData(args.outsilent, False, False, "binary", core.io.silent.SilentFileOptions())

            self.outsilent = args.outsilent

        self.pdb = False
        if not args.pdbdir == '':
            self.pdb = True

            self.pdbdir    = args.pdbdir
            self.outpdbdir = args.outpdbdir

            self.struct_iterator = glob.glob(os.path.join(args.pdbdir, '*.pdb'))

            # Parse the runlist and determine which structures to process
            if args.runlist != '':
                with open(args.runlist, 'r') as f:
                    self.runlist = set([line.strip() for line in f])

                    # Filter the struct iterator to only include those in the runlist
                    self.struct_iterator = [struct for struct in self.struct_iterator if os.path.basename(struct).split('.')[0] in self.runlist]

                    print(f'After filtering by runlist, {len(self.struct_iterator)} structures remain')

        # Assert that either silent or pdb is true, but not both
        assert(self.silent ^ self.pdb), f'Both silent and pdb are set to {args.silent} and {args.pdb} respectively. Only one of these may be active at a time'

        # Setup checkpointing
        self.chkfn = args.checkpoint_name
        self.finished_structs = set()

        if os.path.isfile(self.chkfn):
            with open(self.chkfn, 'r') as f:
                for line in f:
                    self.finished_structs.add(line.strip())

    def record_checkpoint(self, tag):
        '''
        Record the fact that this tag has been processed.
        Write this tag to the list of finished structs
        '''
        with open(self.chkfn, 'a') as f:
            f.write(f'{tag}\n')

    def iterate(self):
        '''
        Iterate over the silent file or pdb directory and run the model on each structure
        '''

        # Iterate over the structs and for each, check that the struct has not already been processed
        for struct in self.struct_iterator:
            tag = os.path.basename(struct).split('.')[0]
            if tag in self.finished_structs:
                print(f'{tag} has already been processed. Skipping')
                continue

            yield struct

    def dump_pose(self, pose, tag):
        '''
        Dump this pose to either a silent file or a pdb file depending on the input arguments
        '''
        if self.pdb:
            # If the outpdbdir does not exist, create it
            # If there are parents in the path that do not exist, create them as well
            if not os.path.exists(self.outpdbdir):
                os.makedirs(self.outpdbdir)

            pdbfile = os.path.join(self.outpdbdir, tag + '.pdb')
            pose.dump_pdb(pdbfile)
        
        if self.silent:
            struct = self.sfd_out.create_SilentStructOP()
            struct.fill_struct(pose, tag)

            self.sfd_out.add_structure(struct)
            self.sfd_out.write_silent_struct(struct, self.outsilent)

    def load_pose(self, tag):
        '''
        Load a pose from either a silent file or a pdb file depending on the input arguments
        '''

        if not self.pdb and not self.silent:
            raise Exception('Neither pdb nor silent is set to True. Cannot load pose')

        if self.pdb:
            pose = pose_from_pdb(tag)
        
        if self.silent:
            pose = Pose()
            self.sfd_in.get_structure(tag).fill_pose(pose)
        
        return pose

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--silent')
    parser.add_argument('--pdbdir')
    parser.add_argument('--out')
    parser.add_argument('--chain', default='B')
    parser.add_argument('--cutoff', default='4.0')
    args = parser.parse_args()


    designs_files = StructManager(args)
    data = ""


    for design in designs_files:
        fixed_res = hotspot_residues(design, args.chain, atom_distance_cutoff=float(args.cutoff))
        print(fixed_res, type(fixed_res))
        fixed_res_string = ''

        for i, resi in enumerate(fixed_res):
            resn = fixed_res[resi]
            fixed_res_string+=f'{resi}{resn},'
        fixed_res_string = fixed_res_string[-1]
        data[design]= fixed_res_string + '\n'

    with open(args.out, 'w+') as f:
        json.dump(data, f)