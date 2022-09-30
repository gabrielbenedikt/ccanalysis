#pragma once
#include <vector>
#include <cstdint>
#include <sstream>
#include <fstream>

#include "histogramset.pb.h"

struct histograms {
    std::vector<std::vector<long long>> offsets;
    std::vector<std::vector<uint32_t>> cc;
    std::vector<std::vector<std::vector<long long>>> cc_tags;
    std::vector<std::vector<uint16_t>> pattern;
    std::vector<double> meastime;
    double resolution;
};

int ser_histogram_protobuf(const histograms* data, const std::string &fname);
histograms deser_histogram_protobuf(const std::string &fn);
