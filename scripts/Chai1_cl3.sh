#!/bin/bash
#SBATCH --job-name=binder_design_BIOC032664
#SBATCH --nodes=1
#SBATCH --gres=gpu:rtx_3090:1
#SBATCH --partition gpu
#SBATCH --mem=18G
#SBATCH --account=BIOC032664

# GPUs A100; node 29, 16GB vRAM, rtx_3090; node 30, a100_sxm4; node 31, 3090 25GB vRAM, V100 32GB vRAM

#SBATCH -e /user/home/vi21227/code/vi21227/code/ProteinDesign/log/Chai1/%A_err.txt
#SBATCH -o /user/home/vi21227/code/vi21227/code/ProteinDesign/log/Chai1/%A_out.txt

# Example running this command - sbatch Chai1.sh -name "test1" -inpdir ../inputs/fc_shaved/  -outdir ../outputs/

# Define Variables


work_dir=/user/home/vi21227/code/vi21227/code/ProteinDesign/outputs/Chai1/$2
input_file=/user/home/vi21227/code/vi21227/code/ProteinDesign/inputs/Chai1/$1/
constraint_file=/user/home/vi21227/code/vi21227/code/ProteinDesign/inputs/Chai1/$1/$3 
chai_script=/user/home/vi21227/code/vi21227/code/ProteinDesign/scripts/Chai1_cyclic_peptide.py

mkdir $work_dir/


echo 'inp' $input_file 'out_dir' $work_dir 'constraint_file' $constraint_file 'chai_script' $chai_script



#Chai1 Hyperparameters

# Contruct Envrioment
module load cuda/12.4.1  
source /user/home/vi21227/.initMamba.sh
conda activate chai
echo "Running Chai1..."
nvcc --version
nvidia-smi


#Chai
echo "Running Chai..."

python $chai_script --fasta_path $input_file \
    --output_dir $work_dir\
    --constraint_path $constraint_file


#Filter designs

