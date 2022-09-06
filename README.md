# several tools to process timetag data

## tsv2hdf 
Convert TSV (tag separated value) files to HDF5.  
HDF5 dataset paths are set in a way to ensure one-way compatibility with how pandas writes HDF5 files (pandas won't be able to read files created by this application. But other applications in this repo are compatible with HDF5 files created by pandas)

## hdf2tsv
The other way round.

## histogram
Reads TSV, HDF5, (zstandard compressed) Cap'n'Proto files and creates coincidence histograms.  
You can configure it by setting values in the file **histogram.cfg**, e.g.   
```
num_threads=8
hist_start=20992
hist_stop=22912
hist_step=6
patterns={{21, 49}, {29,49}}
truncate=0
tagger_resolution=0.015625
``` 
num_threads
: Histogram generation is OpenMP multithreaded. This defines the number of threads to use  

tagger_resolution
: Your timetagger resolution in nanoseconds  

hist_start, hist_stop, hist_step
: Beginning, end, and resolution of histogram (delay between two channels)  

truncate=n
: Only use the first n seconds present in your tagfile for analysis  

patterns
: List of coincidence patterns for which histograms will be created. Each pattern itself is a list of two channel numbers  

# Plot histogram
See scripts/histogram.py for example usage in python.

# How to build

## prerequesties
- HDF5
- Boost
- Protobuf
- CapnProto
- fmtlib

mkdir build  
cd build  
cmake ..  
make  
