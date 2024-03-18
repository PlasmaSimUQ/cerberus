#!/bin/bash

#The file will install some of the commeon missing dependencies, and ignore them if they are already installed. If there
# is a dependency missing after you run this (compilation fails after running this, due to a missing package) then 
# post an issue on the git page please. 

pkgs=(g++ python3 python-is-python3 mpich libboost-all-dev libblas-dev liblapack-dev libmotif-dev libxext-dev libxpm-dev libreadline-dev)
sudo apt-get -y --ignore-missing install "${pkgs[@]}" 

