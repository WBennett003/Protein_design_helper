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

ProteinMPNN_path=$Protein_Design_Helper/tools/LigandMPNN
working_dir=$Protein_Design_Helper/working/SolMPNN/

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


mkdir -p $working_dir/outputs/solMPNN/$2/analysis/

python $working_dir/scripts/analysis/SolMPNN_analysis.py \
    --fasta_filepath "$working_dir/outputs/solMPNN/$2/seqs/$1" \
    --output_dir $working_dir/outputs/solMPNN/$2/analysis/ \
    --chain $8
