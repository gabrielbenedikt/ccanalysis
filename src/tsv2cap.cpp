#include "tsv2cap.h"

namespace bpo = boost::program_options;

int main(int argc, char **argv) {
    
    std::vector<std::string> fnames = {};
    bool compress = false;
    int8_t compression_level = DEFAULT_COMPRESSION_LEVEL;
    /*
     * parse arguments
     */
    try {
        bpo::options_description args("Arguments");
        args.add_options()
            ("help,h", "prints this message")
            ("compress,c", bpo::bool_switch()->default_value(false), "zstd-compress output file")
            ("compression-level,l", bpo::value<int8_t>()->default_value(DEFAULT_COMPRESSION_LEVEL), "zstd compression level (-7...22)")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        bpo::options_description args_hidden("Hidden options");
        args_hidden.add_options()
            ("input-files", bpo::value<std::vector<std::string>>(), "TSV files to convert")
            ;
        
        bpo::options_description cmdline_options;
        cmdline_options.add(args).add(args_hidden);
        
        bpo::positional_options_description args_pos;
        args_pos.add("input-files", -1);
        
        bpo::variables_map args_vm;
        store(bpo::command_line_parser(argc, argv).
              options(cmdline_options).positional(args_pos).run(), args_vm);
        notify(args_vm);
        
        if (args_vm.count("help") != 0) {
            std::cout << "Usage: ./cap2tsv input-files\n";
            std::cout << "Converts TSV files capnproto tag files.\n";
            std::cout << args << "\n";
            return 0;
        }

        if (args_vm.count("input-files") != 0) {
            fnames = args_vm["input-files"].as< std::vector<std::string> >();
        } else {
            std::cout << "no input files specified.\n";
            return 1;
        }
        if (args_vm.count("compression-level") != 0) {
            if (args_vm["compression-level"].as<int8_t>() < MIN_COMPRESSION_LEVEL ) {
                compression_level = MIN_COMPRESSION_LEVEL;
            } else if (args_vm["compression-level"].as<int8_t>() > MAX_COMPRESSION_LEVEL ) {
                compression_level = MAX_COMPRESSION_LEVEL;
            } else {
                compression_level = args_vm["compression-level"].as<int8_t>();
            }
        compress = args_vm["compress"].as<bool>();
        
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
    
    /*
     * process
     */
    for (const auto &fname: fnames) {
        auto starttime = std::chrono::high_resolution_clock::now();
        std::string fn_cap = fname;
        stringreplace(fn_cap, std::filesystem::path(fname).extension(), ".tags");
        std::cout << fname << "\treading..." << std::flush;
        std::vector<long long> data;
        readTSVtags(fname, data);
        std::cout << "ok\twriting..." << std::flush;
        writecapnptags(fn_cap, data, compress, compression_level);
        auto endtime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endtime - starttime).count();
        std::cout << "ok\t" << duration << "ms" << std::endl;
    }
    return 0;
}
