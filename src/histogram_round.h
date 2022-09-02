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
#include "histogramset.pb.h"

/*
 * TODO
 * 
 * don't convert to double
 */

// CONFIG VALUES
struct Config {
    double CS = 0.015625;
    double HIST_START = 0;
    double HIST_STOP  = 0;
    double HIST_STEP  = 0.1;
    int TRUNCATE_S = -1;
    uint16_t NUM_THREADS = 1;
    std::vector<std::vector<uint16_t>> patterns = {};
} cfg;

double ETA = 0;
double SPENT_TIME = 0;

struct histogram_onepattern {
    std::vector<double> offsets;
    std::vector<long> cc;
    std::vector<std::vector<double>> cc_tags;
    std::vector<uint16_t> pattern;
    double meastime;
};

struct histograms {
    std::vector<std::vector<double>> offsets;
    std::vector<std::vector<long>> cc;
    std::vector<std::vector<std::vector<double>>> cc_tags;
    std::vector<std::vector<uint16_t>> pattern;
    std::vector<double> meastime;
};

void separate_tags_per_channels(const long long* tags, const long long numtags, std::vector<std::vector<double>>& tags_per_channel, const std::vector<uint16_t> channels);

// set intersect
histogram_onepattern histogram(const std::vector<std::vector<double>>* tags_per_channel, const std::vector<uint16_t>* channels, const std::vector<double>* offsets, const std::vector<uint16_t> pattern);  

void histograms_to_struct(const std::vector<histogram_onepattern> *pts, histograms *hs);
int histstruct_protobuf_todisk(const histograms* data, const std::string fname);

std::vector<std::string> get_new_tagfiles();
void read_config();


void print_histogram_onepattern(histogram_onepattern h);
#endif // HISTOGRAM_H
