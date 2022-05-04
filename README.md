# RBE

## simulation 
Geant4 simulation of a segment of chromatin fibre in a water volume. Energy deposits and OH reactions are scored for analysis with runClustering. 

run with ./rbe and options as required:
 - -mac filename - to set the input file
 - -out filename - to output a root file
 - -ref - to use the reference geometry for photon beams with 2mm build up, corresponding input file is Co60.in
 - -seed - to set the seed
 - -gui - for visualisation
 - -chemOFF - to run without chemistry

## Clustering 
Analyses energy depositions and OH reactions from the simulation to count strand breaks (SSB, cSSB, DSB).

to run:
- Create/activate a python environment with the same python version as pyroot (root-config --python-version), and install pybind11, numpy and scipy. (Blue pebble use apps/root/6.26.00, gcc 7.5.0, lang/python/miniconda/3.9.7)
- change CMakeLists.txt pybind11_DIR and header directories to the environment.
- build pyClustering.cc 
- change input/output filenames in runClustering.py before running

input.csv should be in the format:

- rootFilename.root,EnergyMeV,LET
