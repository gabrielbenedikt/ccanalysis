#ifndef CCANALYSIS_INTTIME_H
#define CCANALYSIS_INTTIME_H

#include <algorithm>
#include <bits/stdc++.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <H5Cpp.h>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

#include "ccset.pb.h"

#define CS          0.015625        //clockcycle to ns
#define PI          3.141592654
//#define HIST_START  3550.0
//#define HIST_STOP    3590.0
//#define HIST_STEP   0.1

// CONFIG VALUES
long BIN_START  = 0;
long BIN_STOP   = 0;
long BIN_WIDTH  = 7;
int CHN_T = 49;
int CHN_H = 25;
int CHN_V = 29;


double ETA = 0;
double SPENT_TIME = 0;

using namespace std;
using namespace H5;
using namespace std::chrono;

namespace fs = std::filesystem;

struct cc_point {
    long offset = 0;
    int cc_h = 0;
    int cc_v = 0;
    vector<long> cc_h_tags = {};
    vector<long> cc_v_tags = {};
};

struct ccstruct {
    vector<long> offsets;
    vector<long> cc_h;
    vector<long> cc_v;
    vector<vector<long>> cc_h_tags;
    vector<vector<long>> cc_v_tags;
    double meastime;
};

long* readHDFfile(string fn, string dsetpath, long& out_data_len);
void separate_tags_per_channels(long* tags, long numtags, vector<long>& out_countst, vector<long>& out_countsh, vector<long>& out_countsv);
vector<cc_point> find_coincidences(vector<long>* ttags, vector<long>* htags, vector<long>* vtags, vector<long>* offsets);  // stl bst

vector<cc_point> find_coincidences2(vector<long>* ttags, vector<long>* htags, vector<long>* vtags, vector<long>* offsets); // self implemented bst
vector<cc_point> find_coincidences3(vector<long>* ttags, vector<long>* htags, vector<long>* vtags, vector<long>* offsets); // unordered set
vector<cc_point> find_coincidences4(vector<long>* ttags, vector<long>* htags, vector<long>* vtags, vector<long>* offsets); // linear search

void print_ccstruct(ccstruct* dat);
void cc_point_to_ccstruct(vector<cc_point> *pts, ccstruct *ccs);
int ccstruct_protobuf_todisk(ccstruct* data, string fname);
int ccstruct_protobuf_fromdisk(ccstruct* data, string fname);

// helpers
template<typename T>
std::vector<T> arange(T start, T stop, T step = 1);
bool stringreplace(string& str, const string& from, const string& to);

vector<string> get_new_tagfiles();
void read_config();
long roundUp(long numToRound, long multiple);
#endif //CCANALYSIS_INTTIME_H
