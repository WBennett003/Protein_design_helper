#!/bin/bash
#SBATCH --job-name=RfDiff_BIOC032664 
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --partition gpu
#SBATCH --account=BIOC032664

#SBATCH -e /user/home/vi21227/work/code/Projects/Designs/gC1q/DR_like/DR5/logs/%A_err.txt
#SBATCH -o /user/home/vi21227/work/code/Projects/Designs/gC1q/DR_like/DR5/logs/%A_out.txt

#Software paths
interface_scripts_path=/user/home/vi21227/work/code/Projects/Protein_design_helper/tools/dl_binder_design/mpnn_fr/
mpnn_fr_path=/user/home/vi21227/work/code/Projects/Protein_design_helper/tools/dl_binder_design
silent_tools_path=/user/home/vi21227/work/code/Projects/Protein_design_helper/tools/dl_binder_design/mpnn_fr/silent_tools
AF2_inital_guess_path=/user/home/vi21227/work/code/Projects/Protein_design_helper/tools/dl_binder_design/af2_initial_guess
SolubleMPNN_path=/user/home/vi21227/code/vi21227/code/Models/LigandMPNN
chai_script_path=/user/home/vi21227/work/code/Projects/Designs/gC1q/DR_like/DR5
binder_score_path=/user/home/vi21227/work/code/Projects/Designs/gC1q/DR_like/DR5/scripts/

# Define Enviroments Variables
source /user/home/$(whoami)/.bashrc
source /user/home/$(whoami)/.initMamba.sh
echo 'user : '"$(whoami)"

complete_target_seq='' #the comeplete target no shaving
num_struct_designs=10
num_seq_designs=10
omit_AA="C"
binder_chain="A"
interface_cutoff=3.0
solMPNN_temperature=0.5

chai_diff_steps=100
chai_trunk_recycles=3

input_pdb=/user/home/vi21227/work/code/Projects/Designs/gC1q/DR_like/DR5/inps/DR5_AB_1pk6_1DU3_very_shaved.pdb 
output_dir=/user/home/vi21227/work/code/Projects/Designs/gC1q/DR_like/DR5/outs/DR5_p15_100_shaved_withseq
file_prefix='250124_DR5_AB_pt15_100_withseq'

# Contruct Envrioment
echo "Loading RFDiff Env"
module load cuda/12.4.1  
conda activate SE3nv 

nvcc --version
nvidia-smi

#Look at sokrypton/rfdiffusion Github for Documentation - https://github.com/sokrypton/RFdiffusion
# python /user/home/vi21227/work/code/Projects/Protein_design_helper/tools/RFdiffusion/run_inference.py \
#     "inference.input_pdb=$input_pdb" \
#     "inference.output_prefix=$output_dir/RfDiff/$file_prefix" \
#     "contigmap.contigs=[42-42/0 B43-80/0 C81-126]" \
#     'contigmap.provide_seq=[43-125]' \
#     inference.num_designs=$num_struct_designs \
#     "ppi.hotspot_res=[B57,B58,B70,C100,C101]" \
#     diffuser.partial_T=20 \

echo "Loading ProteinMPNN + FastRelax Env"
#Run ProteinMPNN with FastRelax
conda activate proteinmpnn_binder_design

#get fixed labels
echo "adding Fixing Residues labels in PDB"

pdbdir=${output_dir}/RfDiff/
trbdir=${output_dir}/RfDiff/
python $mpnn_fr_path/helper_scripts/addFIXEDlabels.py --pdbdir $pdbdir --trbdir $trbdir --verbose
ls "$trbdir" | grep '\.pdb$' | xargs -I {} rm "$trbdir/{}"

# #convert the generated pdbs onto silents files to save space
echo "Merging PDBs into Silent file"
$silent_tools_path/silentfrompdbs $pdbdir > $output_dir/RfDiff/RfDiff_designs.silent
ls "$pdbdir" | grep '\.trb$' | xargs -I {} rm "$pdbdir/{}"

echo "Running Interface Design"
python $mpnn_fr_path/mpnn_fr/dl_interface_design.py \
    --silent $output_dir/RfDiff/RfDiff_designs.silent \
    --outsilent $output_dir/protein_mpnn_design.silent 

#AF2 interface prediction
echo "AF2 predicting interface"
conda activate af2_binder_design
python $AF2_inital_guess_path/predict.py -silent $output_dir/protein_mpnn_design.silent  -outsilent $output_dir/AF2/af2.silent

#extract interface residues from AF2 prediction
echo "Extracting interface residues"
conda activate BindCraft
python $interface_scripts_path/binder_solmpnn_prep.py --silent $output_dir/af2.silent --out $output_dir/interface_fix.json --chain $binder_chain --cutoff $interface_cutoff

#Run SolubleMPNN on non-interface residues
echo "Redesigning non-interating residues with SolubleMPNN"
conda activate SolMPNN

#Look at LigandMPNN Github for Documentation
python $SolubleMPNN_path/run.py \
    --model_type "soluble_mpnn" \
    --checkpoint_soluble_mpnn "$ProteinMPNN_path/model_params/solublempnn_v_48_002.pt" \
    --seed 42 \
    --pdb_path_multi $output_dir/ \
    --out_folder $output_dir/SolMPNN/ \
    --batch_size $num_seq_designs \
    --number_of_batches 1 \
    --fixed_residues_multi $output_dir/interface_fix.json \
    --omit_AA $OmitAA \
    --chains_to_design $binder_chain \
    --temperature $solMPNN_temperature

echo "Merging Designs"
#Create fasta for chai1 with binder and whole structure
python /user/home/vi21227/work/code/Projects/Designs/gC1q/DR_like/DR5/scripts/fasta_binder_target.py \
    --fasta_file $output_dir/SolMPNN/seqs/ \
    --target_seq $target_seq \
    --chain $binder_chain \
    --output_file $output_dir/SolMPNN/seqs/merged.fasta


#Run Chai predicition
echo "Prediciting final structure with Chai1"
conda activate Chai
python $chai_script_path/Chai1_batch.py \
    --fasta_path $output_dir/SolMPNN/merged.fasta \
    --output_dir $output_dir/Chai/ \
    --num_trunk_recycles $chai_trunk_recycles \
    --num_diffn_timesteps $chai_diff_steps

echo "Getting Design Data"
python $binder_score_path/binder_scoring.py \
    --pdbdir $output_dir/Chai/ \
    --output $output_dir/analysis/binder_stats.csv
