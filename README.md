# several tools to process timetag data

## tsv2hdf 
TSV (tag separated value) files to HDF5

## hdf2tsv
The other way round

## histogram
reads TSV, HDF5, (zstandard compressed) Cap'n'Proto files and creates coincidence histograms
see cc.cfg_example for configuration

## rewinding_ccanalysis
Project specific binary. Won't be usefull for you
Expects HDF5 files as created by the scripts/txttohdf.py script, which itself expects "tab separated values" (TSV) files as created by the UQDevices timetagger.  
Creates protobuf files that are easily readable by python.  

See scripts/histogram.py for example usage in python.

# build
mkdir build  
cd build  
cmake ..  
make  
