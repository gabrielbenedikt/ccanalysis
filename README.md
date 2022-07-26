# coincidence analysis script for rewinding experiment
Expects HDF5 files as created by the scripts/txttohdf.py script, which itself expects "tab separated values" (TSV) files as created by the UQDevices timetagger.  
Creates Cap'n'Proto files that are easily readable by python.  

See scripts/histogram.py for example usage in python.  

# build
mkdir build  
cd build  
cmake ..  
make  
