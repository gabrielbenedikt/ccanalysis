#!/bin/bash
# Instrument binaries, pgo data to /data/pgo, serial make is important to not confuse the pgo generator
env CXXFLAGS='-pipe -O3 -march=native -fomit-frame-pointer -fprofile-dir=~/Nextcloud/projects/rewinding/data/performance/cpp/pgo -fprofile-generate=~/Nextcloud/projects/rewinding/data/performance/cpp/pgo' cmake . -DCMAKE_BUILD_TYPE=Release                   
make -j 1

# Run instrumented program, generate and write pgo data
./runIt

# Use profile data and feed into gcc, correct for threading counter noise, serial make is important to not confuse the pgo generator
env CXXFLAGS='-pipe -O3 -march=native -fomit-frame-pointer -fprofile-dir=~/Nextcloud/projects/rewinding/data/performance/cpp/pgo -fprofile-use=~/Nextcloud/projects/rewinding/data/performance/cpp/pgo -fprofile-correction' cmake . -DCMAKE_BUILD_TYPE=Release
make -j 1
