# import pyrosettacolabsetup; pyrosettacolabsetup.install_pyrosetta()
import pyrosetta; pyrosetta.init()
from pyrosetta.toolbox.mutants import mutate_residue
from pyrosetta import *
from pyrosetta.rosetta.core.pose import Pose

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import argparse

def write_mutation(file_path, mutant_res, position, ddG):
    with open(file_path, 'w+') as f:
        text = f.read()
        text += f"{mutant_res+str(position)},{position},{mutant_res},{ddG}" #name, pos, mutant_res, dgg
        f.write(text)

def binder_residue_scan(pose, binder_chain, offset=0, residue="A", out_file=''):
  rel = pyrosetta.rosetta.protocols.relax.FastRelax()
  score_function = pyrosetta.rosetta.core.scoring.get_score_function()
  rel.set_scorefxn(score_function)

  ddG = []
  RMSD = []
  TM = []

  wt_score = score_function(pose)
  
  num_chains = pose.num_chains()
  for i in range(1, num_chains+1):  # PyRosetta indices are 1-based
    print(i, pose.pdb_info().chain(i), binder_chain, binder_chain == i)
    if i == binder_chain:
      start = pose.chain_begin(i)
      end = pose.chain_end(i)
      binder_length = int(end) - int(start)
      break

  for i in range(binder_length-offset):
    print(f"running residues {residue} at position {i}")
    pose_mutant = pose.clone()
    mutate_residue(pose_mutant,i+start, residue)
    rel.apply(pose_mutant)
    mutant_score = score_function(pose_mutant)
    ddG.append(wt_score-mutant_score)
    print(wt_score-mutant_score)
    write_mutation(out_file, residue, i, float(wt_score-mutant_score))
  return ddG

def get_chain_seq(pose, chain_id):
    return pose.pdb_info().chain(chain_id).sequence()

def plot_deep_scan(data, res, seq, output_dir=''):
    sns.heatmap(data, xticklabels=seq, yticklabels=res)
    plt.savefig(output_dir+'heatmap.png')

def load_pdb(pdb_dir):
    pose = Pose()
    pose_from_file(pose, pdb_dir)
    return pose

def deep_scan(starting_pose, binder_chain_idx, deep=True, output_dir=""):
    data = {}
    if deep:
        residues = ['A']
    
    for res in residues:
        data[res] = binder_residue_scan(pose, binder_chain_idx, residue=res, out_file=output_dir+'mutants.csv')
    

    plot_deep_scan(data, residues, starting_pose.pdb_info().chain(binder_chain_idx).sequence(), output_dir)
    return data

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_file')
    parser.add_argument('--output_dir', default='/user/home/vi21227/code/vi21227/code/ProteinDesign/outputs/deepscan/')
    parser.add_argument('--chain_idx', default=4)
    args = parser.parse_args()

    pose = load_pdb(args.pdb_file)
    print(pose)
    data = deep_scan(pose, int(args.chain_idx), output_dir=args.output_dir)

    df = pd.DataFrame(data)
    df.to_csv(output_dir+'/data.csv')
