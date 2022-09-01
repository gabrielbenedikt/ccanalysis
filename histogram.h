#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <random>
#include <algorithm>
#include <bits/stdc++.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

#include "tools.h"
#include "io.h"
#include "histogramset_int.pb.h"

// CONFIG VALUES
struct Config {
    double CS = 0.015625;                               // tag resolution in ns
    long long HIST_START = 0;                           // start delay for histogram creation. in terms of tag resolution
    long long HIST_STOP  = 0;                           // stop delay for histogram creation. in terms of tag resolution
    long long HIST_STEP  = 1;                           // granularity of histogram. in terms of tag resolution
    int TRUNCATE_S = -1;                                // truncate data after TRUNCATE_S seconds
    uint16_t NUM_THREADS = 1;                           // number of threads the coincidence analysis uses
    uint64_t WND;                                       // coincidence window in terms of tag resolution
    std::vector<std::vector<uint16_t>> patterns = {};   // vector of coincidence patterns, themselves in vector form
} cfg;

struct histogram_onepattern {
    std::vector<long long> offsets;
    std::vector<long> cc;
    std::vector<std::vector<long long>> cc_tags;
    std::vector<uint16_t> pattern;
    double meastime;
};

struct histograms {
    std::vector<std::vector<long long>> offsets;
    std::vector<std::vector<long>> cc;
    std::vector<std::vector<std::vector<long long>>> cc_tags;
    std::vector<std::vector<uint16_t>> pattern;
    std::vector<double> meastime;
};

void separate_tags_per_channels(const long long* tags, const long long numtags, std::vector<std::vector<long long>>& tags_per_channel, const std::vector<uint16_t> channels);

// set intersect
histogram_onepattern histogram(const std::vector<std::vector<long long>>* tags_per_channel, const std::vector<uint16_t>* channels, const std::vector<long long>* offsets, const std::vector<uint16_t> pattern);  

void histograms_to_struct(const std::vector<histogram_onepattern> *pts, histograms *hs);
int histstruct_protobuf_todisk(const histograms* data, const std::string fname);

std::vector<std::string> get_new_tagfiles();
void read_config();


void print_histogram_onepattern(histogram_onepattern h);
#endif // HISTOGRAM_H
