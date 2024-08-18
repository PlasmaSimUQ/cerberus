
# update your system before beginning
sudo apt-get update 

# ensure curl is installed. it allows you to download files directly from urls
sudo apt-get install curl

# go to  /tmp and download then install conda. /tmp is a temporary storage location 
#for files and directories that are created and accessed during system runtime
cd /tmp
curl -O https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh # you can download a different version from the link below
									      # https://repo.anaconda.com/archive/
# install the file 
bash /tmp/Anaconda3-2020.11-Linux-x86_64.sh

# update current session with installed files 
source ~/.bashrc

#========================================== Notes ============================================#
# start the base conda environment
eval "$(/home/kyriakos/anaconda3/bin/conda shell.bash hook)"

#inspect available environments 
#> conda info --envs

#IFF you selected to have conda autostart on terminal start up and you want to change this then execute this command. 
#> conda config --set auto_activate_base false

#clone environments... i.e. make a copy to some stuff, change the name cloned_name and select the starting 
#environment base other 
# conda create --name cloned_name --clone base
#	Source:      /home/kyriakos/anaconda3
#	Destination: /home/kyriakos/anaconda3/envs/cloned

#Get a yaml file for the environment 
#> conda env export > environment.yml

# Create from yaml file 
#>conda env create -f environment.yml

# create an environment from specific packages ...
conda create -n cerberus_python python=3.10.12 numpy=1.26.4 scipy=1.12.0 matplotlib=3.5.1

conda activate cerberus_python
# ============================================ Misc =============================================#
# check the version and see if it even works
#>conda -v

#update 
#>conda update conda
#>conda update anaconda

#conda create --name py3 python=3

#conda activate py3

#conda deactivate


