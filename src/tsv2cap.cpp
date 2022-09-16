#include "tsv2cap.h"

namespace bpo = boost::program_options;

int main(int argc, char **argv) {
    
    std::vector<std::string> fnames = {};
    bool compress = false;
    uint32_t compression_level = 3;
    /*
     * parse arguments
     */
    try {
        bpo::options_description args("Arguments");
        args.add_options()
            ("help,h", "prints this message")
            ("compress,c", bpo::bool_switch()->default_value(false), "zstd-compress output file")
            ("compression-level,l", bpo::value<uint32_t>()->default_value(3), "zstd compression level")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        bpo::options_description hidden("Hidden options");
        hidden.add_options()
            ("input-files", bpo::value<std::vector<std::string>>(), "TSV files to convert")
            ;
        
        bpo::options_description cmdline_options;
        cmdline_options.add(args).add(hidden);
        
        bpo::positional_options_description p;
        p.add("input-files", -1);
        
        bpo::variables_map vm;
        store(bpo::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);
        
        if (vm.count("help")) {
            std::cout << "Usage: ./cap2tsv input-files\n";
            std::cout << "Converts TSV files capnproto tag files.\n";
            std::cout << args << "\n";
            return 0;
        }

        if (vm.count("input-files")) {
            fnames = vm["input-files"].as< std::vector<std::string> >();
        } else {
            std::cout << "no input files specified.\n";
            return 1;
        }
        if (vm.count("compression-level")) {
            compression_level = vm["compression-level"].as< uint32_t >();
        } else {
            std::cout << "no input files specified.\n";
            return 1;
        }
        compress = vm["compress"].as<bool>();
        
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
    for (auto fn: fnames) {
        auto starttime = std::chrono::high_resolution_clock::now();
        std::string fn_cap = fn;
        stringreplace(fn_cap, std::filesystem::path(fn).extension(), ".tags");
        std::cout << fn << "\treading..." << std::flush;
        long long tsv_len = 0;
        std::vector<long long> data;
        readTSVtags(fn, data, tsv_len);
        std::cout << "ok\twriting..." << std::flush;
        writecapnptags(fn_cap, data, compress, compression_level);
        auto endtime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endtime - starttime).count();
        std::cout << "ok\t" << duration << "ms" << std::endl;
    }
    return 0;
}
