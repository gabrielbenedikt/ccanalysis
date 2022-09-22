#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <filesystem>
#include <boost/program_options.hpp>
#include <chrono>
#include "io.h"

void split_hdf5(const std::string &fn, const uint16_t numparts = 0, const uint32_t numtags = 0);
void split_tsv(const std::string &fn, const uint16_t numparts = 0, const uint32_t numtags = 0);
void split_cap(const std::string &fn, const uint16_t numparts = 0, const uint32_t numtags = 0);
