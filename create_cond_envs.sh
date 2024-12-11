
source ~/.initMaba.sh

# SetUP Chai vitrual env
conda create -n Chai

conda activate Chai

pip install wandb chai-lab

# SetUP SolubleMPNN vitrual env
conda create -n SolMPNN

conda activate Chai

pip install wandb chai-lab

# SetUP RfDiffusion vitrual env
conda create -n RfDiff

conda activate Chai

pip install wandb chai-lab


# SetUP Bindcraft vitrual env
conda create -n BindCraft

