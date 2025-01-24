from pathlib import Path

import pandas as pd
import numpy as np
import torch

from chai_lab.chai1 import run_inference
import os

import argparse

# We use fasta-like format for inputs.
# - each entity encodes protein, ligand, RNA or DNA
# - each entity is labeled with unique name;
# - ligands are encoded with SMILES; modified residues encoded like AAA(SEP)AAA

# Example given below, just modify it

parser = argparse.ArgumentParser(description='Script to run BindCraft binder design.')
parser.add_argument('--fasta_path')

parser.add_argument('--output_dir')
parser.add_argument('--num_trunk_recycles', default=3)
parser.add_argument('--num_diffn_timesteps', default=200)
parser.add_argument('--seed', default=42)
parser.add_argument('--use_esm_embeddings', default=True)

args = parser.parse_args()

fasta_path = args.fasta_path
output_dir = args.output_dir
num_trunk_recycles = args.num_trunk_recycles
num_diffn_timesteps = args.num_diffn_timesteps
seed = args.seed
use_esm_embeddings = args.use_esm_embeddings

print(fasta_path)
print(output_dir)
fastas_files = [p for p in Path(fasta_path).glob('*.fasta')]
print(fastas_files)
# print([x for x in fastas_files])

data = {
    "fasta_name" : [],
    "pAE" : [],
    "pLDDT" : [],
    "pTM" : [],
    "pDE" : [],
    "i_pTM" : [],
}

for fasta in fastas_files:
    with open(fasta, 'r') as f:
        input_text = f.read().strip()
    
    with open(fasta, 'w') as f:
        f.write(input_text)

    f = Path(fasta)
    print(output_dir, str(output_dir)+'/'+str(fasta).split('/')[-1])

    candidates = run_inference(
        fasta_file=f,
        output_dir=Path(str(output_dir)+'/'+str(fasta).split('/')[-1]),
        # 'default' setup
        num_trunk_recycles=num_trunk_recycles,
        num_diffn_timesteps=num_diffn_timesteps,
        seed=seed,
        device=torch.device("cuda:0"),
        use_esm_embeddings=use_esm_embeddings,
    )

    cif_paths = candidates.cif_paths
    scores = [rd.aggregate_score for rd in candidates.ranking_data]


    data['fasta_name'].append(str(fasta).split('/')[-1])
    data['pAE'].append(mean(candidates.pae))
    data['pDE'].append(mean(candidates.pde))
    data['pLDDT'].append(mean(candidates.plddt))
    data['pTM'].append(sum([x.ptm_scores.complex_ptm for x in candidates.ranking_data]))
    data['i_pTM'].append(sum([x.ptm_scores.interface_ptm for x in candidates.ranking_data]))

## Load scores and create Dataframe...

df = pd.DataFrame(data)
df.to_csv(output_dir+'/'+'summary.csv')


