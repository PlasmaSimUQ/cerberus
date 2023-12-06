# =======================================================================
#import matplotlib.pyplot as mpl
import numpy as np 
import os, sys, gc, numpy as np # standard modules
import pylab as plt, matplotlib.gridspec as gridspec, matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from multiprocessing import Pool
import pdb, math
import copy 
import re

visulaisation_code_folder     ="/home/kyriakos/Documents_ubuntu/cerberus/vis"
derived_functions_code_folder ="/home/kyriakos/Documents_ubuntu/cerberus/vis"

if visulaisation_code_folder not in sys.path:
  sys.path.insert(0, visulaisation_code_folder)
if derived_functions_code_folder not in sys.path:
  sys.path.insert(0, derived_functions_code_folder)

import PHM_MFP_Solver_Post_functions_v6 as phmmfp # running version 2 
from get_boxlib import ReadBoxLib, get_files, parse_header

def atoi(text):
  return int(text) if text.isdigit() else text

def natural_keys(text):
  '''
  alist.sort(key=natural_keys) sorts in human order
  http://nedbatchelder.com/blog/200712/human_sorting.html
  (See Toothy's implementation in the comments)
  '''
  return [ atoi(c) for c in re.split('plt', text) ]

file_dir = os.getcwd()

#file_dir ="/home/kyriakos/Documents/Code/000_cerberus_dev/cerberus/Exec/testing/Collisions/GhoshThermalEquilibration_collisions/"

#plt_file = "GhoshThermalEquilibration"

#basepath = file_dir
#for entry in os.listdir(basepath):
#    if os.path.isfile(os.path.join(basepath, entry)):
#        print(entry)

# get a list of all the files in this directory
files = get_files(file_dir, include=["plt"], exclude=["chk", ".png", "inputs"], times=[], tol=1e-4, get_all=True)

if len(files) == 0: pdb.set_trace()

files.sort(key=natural_keys)

visit_file_object = open(file_dir + "/DataDirectory.visit", "w")

for directoryName in files:
  visit_file_object.write(directoryName + "/Header\n")

visit_file_object.close()




