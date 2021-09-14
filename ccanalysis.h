#ifndef CCANALYSIS_H
#define CCANALYSIS_H

#define CC_ALGO STL_SETINTERSECT
                //STL_SETINTERSECT_TRIGGER
                //BST
                //UNORDERED_SET
                //LINEAR_SEARCH
                //STL_BST

#include <algorithm>
#include <bits/stdc++.h>
#include <chrono>
//#include <execution>
#include <filesystem>
#include <fstream>
#include <H5Cpp.h>
#include <iterator>
#include <iostream>
#include <omp.h>
#include <string>
// #include <unordered_set>
#include <vector>

#include "ccset.pb.h"

#define CS          0.015625        //clockcycle to ns
#define PI          3.141592654

// CONFIG VALUES
double HIST_START = 0;
double HIST_STOP  = 0;
double HIST_STEP  = 0.1;
double FPGA_DELAY_MIN = 0;
double FPGA_DELAY_MAX = 0;
int CHN_TR = 49;
int CHN_FPGA = 1;
int CHN_H = 25;
int CHN_V = 29;
int TRUNCATE_S = -1;

double ETA = 0;
double SPENT_TIME = 0;

struct cc_point {
    double offset = 0;
    int cc_h = 0;
    int cc_v = 0;
    std::vector<double> cc_h_tags = {};
    std::vector<double> cc_v_tags = {};
};

struct ccstruct {
    std::vector<double> offsets;
    std::vector<long> cc_h;
    std::vector<long> cc_v;
    std::vector<std::vector<double>> cc_h_tags;
    std::vector<std::vector<double>> cc_v_tags;
    long num_trigger_tags = 0;
    long num_fpga_tags = 0;
    double meastime;
};

long long* readHDFfile(const std::string fn, const std::string dsetpath, long long& out_data_len);
void separate_tags_per_channels(const long long* tags, const long long numtags, std::vector<double>& out_countst, std::vector<double>& out_countsh, std::vector<double>& out_countsv, std::vector<double>& out_countsfpga);

#if CCALGO==STL_SETINTERSECT
// set intersect
std::vector<cc_point> find_coincidences(const std::vector<double>* ttags, const std::vector<double>* htags, const std::vector<double>* vtags, const std::vector<double>* offsets);  
#elif CCALGO==STL_SETINTERSECT_TRIGGER
//set intersect. keep only coincidences where fpga is present too
std::vector<cc_point> find_coincidences(const std::vector<double>* ttags, const std::vector<double>* htags, const std::vector<double>* vtags, const std::vector<double>* offsets, const std::vector<double>* ftags); 
#elif CCALGO==BST
// self implemented bst
std::vector<cc_point> find_coincidences(const std::vector<double>* ttags, const std::vector<double>* htags, const std::vector<double>* vtags, const std::vector<double>* offsets); 
#elif CCALGO==UNORDERED_SET
// unordered set
std::vector<cc_point> find_coincidences(const std::vector<double>* ttags, const std::vector<double>* htags, const std::vector<double>* vtags, const std::vector<double>* offsets); 
#elif CCALGO==LINEAR_SEARCH
// linear search
std::vector<cc_point> find_coincidences(const std::vector<double>* ttags, const std::vector<double>* htags, const std::vector<double>* vtags, const std::vector<double>* offsets);
#elif CCALGO==STL_BST
// stl bst
std::vector<cc_point> find_coincidences(const std::vector<double>* ttags, const std::vector<double>* htags, const std::vector<double>* vtags, const std::vector<double>* offsets); 
#endif

void print_ccstruct(const ccstruct* dat);
void cc_point_to_ccstruct(const std::vector<cc_point> *pts, ccstruct *ccs);
int ccstruct_protobuf_todisk(const ccstruct* data, const std::string fname);
int ccstruct_protobuf_fromdisk(ccstruct* data, const std::string fname);

// helpers
template<typename T>
std::vector<T> arange(const T start, const T stop, const T step = 1);
bool stringreplace(std::string& str, const std::string& from, const std::string& to);

std::vector<std::string> get_new_tagfiles();
void read_config();

#endif //CCANALYSIS_H
