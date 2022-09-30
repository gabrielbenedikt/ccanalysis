#include "histcat.h"
#include "io_protobuf.h"

namespace bpo = boost::program_options;

enum APPEND_MODE {
    PATTERN = 1,
    OFFSET = 2,
    CC = 3
};

int main(int argc, char **argv)
{ 
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    int8_t mode = -1;
    std::string modestr;
    std::string outfn;
    /*
     * get files from prompt
     */
    std::vector<std::string> histfiles;
    try {
        bpo::options_description args("Arguments");
        args.add_options()
            ("help,h", "prints this message")
            ("input-files", bpo::value<std::vector<std::string>>(), "Histogram output files to combine")
            ("append,a", bpo::value<std::string>()->default_value("tag"), "Datatype to append. <tag> or <pattern> or <offset>.")
            ("out,o", bpo::value<std::string>(), "Output file name.")
            ;
        bpo::options_description cmdline_options;
        cmdline_options.add(args);
        bpo::positional_options_description args_pos;
        args_pos.add("input-files", -1);
        
        bpo::variables_map args_vm;
        store(bpo::command_line_parser(argc, argv).
              options(cmdline_options).positional(args_pos).run(), args_vm);
        notify(args_vm);
        
        if (args_vm.count("help") != 0) {
            std::cout << "Usage: ./histcat -a=[tag | pattern | offset] input-files\n";
            std::cout << "Combines histogram output files.\n";
            std::cout << args << "\n";
            return 0;
        }

        if (args_vm.count("input-files") != 0) {
            histfiles = args_vm["input-files"].as< std::vector<std::string> >();
        } else {
            std::cout << "no input files specified." << std::endl;
            return 2;
        }
        if (histfiles.size() < 2) {
            std::cerr << "Need more than 1 file to combine them..." << std::endl;
        }
        
        if (args_vm.count("out") != 0) {
            outfn = args_vm["out"].as< std::string >();
        } else {
            auto p = std::filesystem::path(histfiles[0]);
            outfn = p.stem().string() + "_combined" + p.extension().string();
        }
        
        if (args_vm.count("append") != 0) {
            modestr = toLower(args_vm["append"].as< std::string >());
        } else {
            std::cout << "No mode of operation specified." << std::endl;
            return 2;
        }
        
        if (modestr == "pat" || modestr == "pattern") {
            mode = APPEND_MODE::PATTERN;
        } else if (modestr == "os" || modestr == "offset") {
            mode = APPEND_MODE::OFFSET;
        } else if (modestr == "cc" || modestr == "tag") {
            mode = APPEND_MODE::CC;
        } else {
            std::cerr << "append mode must be in [pat | pattern, os | offset, cc | tags]" << std::endl;
            return 3;
        }
        
        
    }
    catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
        return 1;
    }
    
    std::vector<histograms> vhistograms;
    for (const auto &fn : histfiles) {
        histograms h = deser_histogram_protobuf(fn);
        vhistograms.emplace_back(h);
    }
    
    histograms h_combined;
    switch (mode) {
        case APPEND_MODE::PATTERN:
            histograms_append_pattern(h_combined, vhistograms);
            ser_histogram_protobuf(&h_combined, outfn);
            break;
        case APPEND_MODE::OFFSET:
            if (histogram_append_offsets(h_combined, vhistograms) == 0) {
                ser_histogram_protobuf(&h_combined, outfn);
            } else {
                std::cerr << "Appending tags failed." << std::endl;
            }
            break;
        case APPEND_MODE::CC:
            if (histogram_append_cctags(h_combined, vhistograms) == 0) {
                ser_histogram_protobuf(&h_combined, outfn);
            } else {
                std::cerr << "Appending tags failed." << std::endl;
            }
            break;
        default:
            std::cerr << "Undefined append mode " << mode << std::endl;
            break;
    }
    
    google::protobuf::ShutdownProtobufLibrary();
    
    return 0;
}

void histograms_append_pattern(histograms &hout, std::vector<histograms> &hvec) {
    hout.resolution = hvec[0].resolution;
    for (const auto &h: hvec) {
        hout.meastime.insert(hout.meastime.end(), h.meastime.begin(), h.meastime.end());
        hout.pattern.insert(hout.pattern.end(), h.pattern.begin(), h.pattern.end());
        hout.cc.insert(hout.cc.end(), h.cc.begin(), h.cc.end());
        hout.offsets.insert(hout.offsets.end(), h.offsets.begin(), h.offsets.end());
        hout.cc_tags.insert(hout.cc_tags.end(), h.cc_tags.begin(), h.cc_tags.end());
    }
}

int histogram_append_cctags(histograms &hout, std::vector<histograms> &hvec) {
    hout = hvec[0];
    //for all remaining structs
    for (size_t i=1; i<hvec.size(); ++i) {
        //check if pattern matches
        bool found = false;
        for (size_t j=0; j<hout.pattern.size(); ++j) {
            for (size_t k=0; k<hvec[i].pattern.size(); ++k) {
                if (hout.pattern[j] == hvec[i].pattern[k]) {
                    found = true;
                    hout.meastime[j] += hvec[i].meastime[k];
                    
                    //check if offsets are the same
                    if (hout.offsets[j] == hvec[i].offsets[k]) {
                        for (size_t os = 0; os < hvec[i].offsets[k].size(); ++os) {
                            hout.cc_tags[j][os].insert(hout.cc_tags[j][os].end(),
                                                                hvec[i].cc_tags[k][os].begin(),
                                                                hvec[i].cc_tags[k][os].end());
                            hout.cc[j][os] += hvec[i].cc[k][os];
                        }
                    } else {
                        std::cerr << "offset mismatch " << std::endl;
                        return -1;
                    }
                }
            }
        }
        if (!found) {
            std::cerr << "Error: could not find same pattern in files" << std::endl;
            return -1;
        }
    }
    
    for (auto &tvv: hout.cc_tags) {
        for (auto &tv: tvv) {
            std::sort(tv.begin(), tv.end());
        }
    }
    
    return 0;
}

int histogram_append_offsets(histograms &hout, std::vector<histograms> &hvec) {
    hout = hvec[0];
    for (size_t i = 1; i<hvec.size(); ++i) {
        //check if pattern matches
        bool found = false;
        for (size_t j=0; j<hout.pattern.size(); ++j) {
            for (size_t k=0; k<hvec[i].pattern.size(); ++k) {
                if (hout.pattern[j] == hvec[i].pattern[k]) {
                    found = true;
                    
                    for (size_t osi=0; osi<hvec[i].offsets[k].size(); ++osi) {
                        hout.offsets[j].emplace_back(hvec[i].offsets[k][osi]);
                        hout.cc[j].emplace_back(hvec[i].cc[k][osi]);
                        hout.cc_tags[j].emplace_back(hvec[i].cc_tags[k][osi]);
                    }
                }
            }
        }
        if (!found) {
            std::cerr << "Error: could not find same pattern in files" << std::endl;
            return -1;
        }
    }
    
    //sort
    for (size_t i = 0; i < hout.pattern.size(); ++i) {
        std::vector<size_t> idx = get_sort_perm(hout.offsets[i]);
        hout.offsets[i] = vector_permute(hout.offsets[i],idx);
        hout.cc[i] = vector_permute(hout.cc[i],idx);
        hout.cc_tags[i] = vector_permute(hout.cc_tags[i],idx);
    }
    return 0;
}
