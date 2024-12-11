#!/bin/bash
#SBATCH --job-name=binder_design_BIOC032664
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --partition gpu
#SBATCH --mem=18G
#SBATCH --account=BIOC032664

#SBATCH -e /user/home/vi21227/code/vi21227/code/ProteinDesign/log/SolubleMPNN/%A_err.txt
#SBATCH -o /user/home/vi21227/code/vi21227/code/ProteinDesign/log/SolubleMPNN/%A_out.txt

module load cuda/12.4.1  
source /user/home/vi21227/.initMamba.sh
conda activate ligandmpnn

ProteinMPNN_path=/user/home/vi21227/code/vi21227/code/Models/LigandMPNN
working_dir=/user/home/vi21227/code/vi21227/code/ProteinDesign

# input_file="cSD1.pdb"
# out_dir="cSD1_run1"

#sbatch -job_name 'test' -inp_file 'MHCI_1_test.pdb' -batch_size 10 -number_of_batches 10 -fix_residues 'A8, A102, B192' -OmitAA 'FAMILYVWC' -chain_to_redesign "E"
echo 'job_name : ' $2 ', inp_file : ' $1 ', batch_size : ' $3 ', number_of_batches : ' $4 ', fix_residues :' $5 ', OmitAA : ' $6 ', Chain_to_redesign : ' $7', chain index' $8'Redesign temperature' $9

python $ProteinMPNN_path/run.py \
    --model_type "soluble_mpnn" \
    --checkpoint_soluble_mpnn "$ProteinMPNN_path/model_params/solublempnn_v_48_002.pt" \
    --seed 42 \
    --pdb_path "$working_dir/inputs/$1" \
    --out_folder "$working_dir/outputs/solMPNN/$2/" \
    --batch_size $3 \
    --number_of_batches $4 \
    --fixed_residues "$5" \
    --omit_AA "$6" \
    --chains_to_design "$7" \
    --temperature $9

# MHCI fixed residues "E5 E7 E24 E25 E26 E59 E62 E63 E64 E65 E66 E67 E68 E69 E70 E71 E72 E73 E74 E75 E76 E77 E78 E79 E80 E81 E85 E87 E95 E96 E97 E98 E99 E114 E116 E117 E118 E123 E124 E130 E140 E142 E143 E144 E146 E147 E148 E149 E150 E151 E152 E153 E154 E155 E156 E157 E158 E159 E160 E161 E162 E163 E163 E164 E166 E167 E168" \
# MHCI hydrophobic bits "E8 E10 E12 E23 E27 E85 E113"
# to get interacting residues in pymol
# OmitAA "FTVIYC"

mkdir -p $working_dir/outputs/solMPNN/$2/analysis/

python $working_dir/scripts/analysis/SolMPNN_analysis.py \
    --fasta_filepath "$working_dir/outputs/solMPNN/$2/seqs/$1" \
    --output_dir $working_dir/outputs/solMPNN/$2/analysis/ \
    --chain $8
