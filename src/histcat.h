#pragma once

#include "io.h"
#include "tools.h"
#include "io_protobuf.h"
#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <iostream>

void histograms_append_pattern(histograms &hout, std::vector<histograms> &hvec);
int histogram_append_cctags(histograms &hout, std::vector<histograms> &hvec);
int histogram_append_offsets(histograms &hout, std::vector<histograms> &hvec);
