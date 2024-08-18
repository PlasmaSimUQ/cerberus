#move this file to e.g. Waveguide and test. You will need to add the first two lines to the run script you want to run or execute them before you run the file of interest 
eval "$(/home/kyriakos/anaconda3/bin/conda shell.bash hook)"

conda activate cerberus_python

python3 plot.py
#python3 movie.py
