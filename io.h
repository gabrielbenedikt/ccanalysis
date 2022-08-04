#ifndef IO_H
#define IO_H

#include <algorithm>
#include <iostream>
#include <string>
#include <string.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <exception>
#include <cassert>
#include <H5Cpp.h>
#include <sys/stat.h>
#include "hdf5.h"
#include "tools.h"

enum HDF5_COMP_ALG {
    HDF5_COMP_ALG_ZLIB = 0,
    HDF5_COMP_ALG_SZIP = 1,
    HDF5_COMP_ALG_NBIT = 2,
    HDF5_COMP_ALG_NONE = 3
};

bool fileExists(const std::string& fn);

void readTSVtags(const std::string fn, std::vector<long long> &result, long long& out_data_len);
long long* readHDF5tags(const std::string fn, long long& out_data_len);

void lltoTSV(const std::string fn, const long long* data, const long long len);
void writeHDFtags(const std::string fn, const std::vector<long long> &r, const uint8_t compression_alg=HDF5_COMP_ALG_ZLIB, const uint8_t compression_level=5);
void writeHDFtagsC(const std::string fn, const std::vector<long long> &r, const uint8_t compression_alg=HDF5_COMP_ALG_ZLIB, const uint8_t compression_level=5);
std::vector<std::string> get_new_tagfiles(const std::string rawext, const std::string analyzed_ext="", const std::vector<std::string> excludes = {});

#endif // IO_H
