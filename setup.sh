
export$PATH=$PATH:$pwd

bash get_model_params.sh "./model_params"

bash create_cond_envs.sh

cd tools

git clone https://github.com/dauparas/LigandMPNN
git clone https://github.com/sokrypton/RFdiffusion 
git clone https://github.com/martinpacesa/BindCraft


cd LigandMPNN
