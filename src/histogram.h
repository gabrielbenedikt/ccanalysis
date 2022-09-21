#pragma once

#include <random>
#include <algorithm>
#include <bits/stdc++.h>
#include <boost/program_options.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>
#include <memory>

#include "tools.h"
#include "io.h"
#include "histogramset.pb.h"

#define SEC_TO_NS 1000000000
#define MIN_TO_SEC 60

// CONFIG VALUES
struct Config {
    double CS            = 0.015625;                    // tag resolution in ns
    long long HIST_START = 0;                           // start delay for histogram creation. in terms of tag resolution
    long long HIST_STOP  = 0;                           // stop delay for histogram creation. in terms of tag resolution
    long long HIST_STEP  = 1;                           // granularity of histogram. in terms of tag resolution
    uint64_t TRUNCATE_S  = 0;                           // truncate data after TRUNCATE_S seconds
    uint16_t NUM_THREADS = 1;                           // number of threads the coincidence analysis uses
    std::vector<std::vector<uint16_t>> patterns = {};   // vector of coincidence patterns, themselves in vector form
};

struct histogram_onepattern {
    std::vector<long long> offsets;
    std::vector<uint32_t> cc;
    std::vector<std::vector<long long>> cc_tags;
    std::vector<uint16_t> pattern;
    double meastime;
};

struct histograms {
    std::vector<std::vector<long long>> offsets;
    std::vector<std::vector<uint32_t>> cc;
    std::vector<std::vector<std::vector<long long>>> cc_tags;
    std::vector<std::vector<uint16_t>> pattern;
    std::vector<double> meastime;
    double resolution;
};

void separate_tags_per_channels(const std::vector<long long> &tags, std::vector<std::vector<long long>>& tags_per_channel, const std::vector<uint16_t> &channels, const Config& cfg);

histogram_onepattern histogram(const std::vector<std::vector<long long>> &tags_per_channel, const std::vector<uint16_t> &channels, const std::vector<long long> &offsets, const std::vector<uint16_t> &pattern, const long long &wnd);

void histograms_to_struct(const std::vector<histogram_onepattern> *pts, histograms *hs, Config& cfg);
int histstruct_protobuf_todisk(const histograms* data, const std::string &fname);

std::vector<std::string> get_new_tagfiles(const Config& cfg);
Config read_config();

void print_histogram_onepattern(histogram_onepattern h);

