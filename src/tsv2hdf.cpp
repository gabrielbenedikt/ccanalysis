#include "tsv2hdf.h"

namespace bpo = boost::program_options;

int main(int argc, char **argv) {

    uint8_t compression_alg = HDF5_COMP_ALG_ZLIB;
    uint8_t compression_level = 6;
    std::vector<std::string> fnames = {};
    /*
     * parse arguments
     */
    
    try {
        bpo::options_description args("Arguments");
        args.add_options()
            ("help,h", "prints this message")
            ("compression-level,L", bpo::value<int>()->default_value(6), "set HDF5 compression level if applicable (0...9). Default 6.")
            ("compression-algorithm,A", bpo::value<std::string>()->default_value("zlib"), "set HDF5 compression algorithm (zlib, szip, nbit, none). Default zlib")
            ;
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
            std::cout << "Usage: ./tsv2hdf [-L int] [-A str] input-files\n";
            std::cout << "Converts tab separated values (TSV) files to HDF5.\n";
            std::cout << args << "\n";
            return 0;
        }

        if (vm.count("compression-level")) {
            if (vm["compression-level"].as<int>() < 0 ) {
                compression_level = 0;
            } else if (vm["compression-level"].as<int>() > 9 ) {
                compression_level = 9;
            } else {
                compression_level = (uint8_t)vm["compression-level"].as<int>();
            }
        }
        
        if (vm.count("compression-algorithm")) {
            std::string algstr = vm["compression-algorithm"].as<std::string>();
            std::transform(algstr.begin(), algstr.end(), algstr.begin(), ::tolower);
            if (algstr == "szip") {
                compression_alg = HDF5_COMP_ALG_SZIP;
            } else if (algstr == "nbit") {
                compression_alg = HDF5_COMP_ALG_NBIT;
            } else if (algstr == "none") {
                compression_alg = HDF5_COMP_ALG_NONE;
            } else {
                compression_alg = HDF5_COMP_ALG_ZLIB;
            }
        }
        
        if (vm.count("input-files")) {
            fnames = vm["input-files"].as< std::vector<std::string> >();
        } else {
            std::cout << "no input files specified.\n";
            return 1;
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
    for (auto fn: fnames) {
        auto starttime = std::chrono::high_resolution_clock::now();
        std::string fn_hdf = fn;
        stringreplace(fn_hdf, std::filesystem::path(fn).extension(), ".h5");
        std::cout << fn << "\treading..." << std::flush;
        std::vector<long long> r;
        readTSVtags(fn, r);
        std::cout << "ok\twriting..." << std::flush;
        writeHDFtags(fn_hdf, r, compression_alg, compression_level);
        auto endtime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endtime - starttime).count();
        std::cout << "ok\t" << duration << "ms" << std::endl;
    }
    return 0;
}
