
source ~/.initMaba.sh

# SetUP Chai vitrual env
conda create -n Chai python=3.11

conda activate Chai

pip install wandb git+https://github.com/chaidiscovery/chai-lab.git

# SetUP SolubleMPNN vitrual env
conda create -n SolMPNN python=3.11

conda activate SolMPNN 

pip install -r $Protein_Design_Helper/tools/LigandMPNN/requirements.txt

pip install wandb 

# SetUP RfDiffusion vitrual env
conda create -n RfDiff python=3.11

conda activate RfDiff 

pip install RfDiff

# SetUP Bindcraft vitrual env
conda create -n BindCraft python=3.11

cd $Protein_Design_Helper/tools/BindCraft
bash install_bindcraft.sh --cuda '12.4' --pkg_manager 'conda'
pip install wandb

#Set up general enviroment

conda create -n bio python=3.11

pip install biopython wandb numpy pandas matplotlib seaborn scikit-learn 
pip install pyrosetta-installer 
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()' #may need key for this from baker