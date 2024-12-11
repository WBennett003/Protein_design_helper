#!/bin/bash
#SBATCH --job-name=binder_design_BIOC032664
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --partition gpu
#SBATCH --mem=18G
#SBATCH --account=BIOC032664

#SBATCH -e /user/home/vi21227/code/vi21227/code/ProteinDesign/log/RfDiffusion/%A_err.txt
#SBATCH -o /user/home/vi21227/code/vi21227/code/ProteinDesign/log/RfDiffusion/%A_out.txt


# Define Variables

input_dir="/user/home/vi21227/work/vi21227/code/ProteinDesign/inputs"
output_path="/user/home/vi21227/code/vi21227/code/ProteinDesign/outputs/RfDiff"
echo 'job name : ' $1 'Input_file : ' $input_dir/$2 'Output_path : ' $output_path/$3 'Num_design : ' $4 'T steps : ' $5 'Contigs :' $6 'partial diffusion steps' $7  

job_name=$1
input_dir=$input_dir/$2
output_path=$output_path/$3
num_designs=$4
iterations=$5
p_iterations=$7
# partial_T=$7
contigs='contigmap.contigs='$6
mkdir $output_path
mkdir $output_path/$1

RfDiffusion_dir=/user/home/vi21227/work/vi21227/code/Models/RFdiffusion

# Contruct Envrioment
module load cuda/12.4.1  
source /user/home/vi21227/.initMamba.sh
conda activate RfDiff 

nvcc --version
nvidia-smi

cd /user/home/vi21227/code/vi21227/code/ProteinDesign/outputs/RfDiff


python $RfDiffusion_dir/run_inference.py $contigs \
    inference.input_pdb=$input_dir \
    inference.output_prefix="$output_path/$job_name" \
    inference.num_designs=$num_designs \
    diffuser.T=$iterations \
    diffuser.partial_T=$p_iterations \
    $8 \
    'ppi.hotspot_res=[A8,A10,A12,A23,A27,A85,A113]'  \
    #inference.symmetry='C2'
    # diffuser.partial_T=10 
#     # inference.ckpt_override_path=models/Complex_beta_ckpt.pt #Beta model


