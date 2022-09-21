#include "tsv2hdf.h"

namespace bpo = boost::program_options;

int main(int argc, char **argv) {

    uint8_t compression_alg = HDF5_COMP_ALG_ZLIB;
    uint8_t compression_level = DEFAULT_COMPRESSION_LEVEL;
    std::vector<std::string> fnames = {};
    /*
     * parse arguments
     */
    
    try {
        bpo::options_description args("Arguments");
        args.add_options()
            ("help,h", "prints this message")
            ("compression-level,L", bpo::value<int>()->default_value(DEFAULT_COMPRESSION_LEVEL), "set HDF5 compression level if applicable (0...9). Default 6.")
            ("compression-algorithm,A", bpo::value<std::string>()->default_value("zlib"), "set HDF5 compression algorithm (zlib, szip, nbit, none). Default zlib")
            ;
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
            std::cout << "Usage: ./tsv2hdf [-L int] [-A str] input-files\n";
            std::cout << "Converts tab separated values (TSV) files to HDF5.\n";
            std::cout << args << "\n";
            return 0;
        }

        if (args_vm.count("compression-level") != 0) {
            if (args_vm["compression-level"].as<int>() < MIN_COMPRESSION_LEVEL ) {
                compression_level = MIN_COMPRESSION_LEVEL;
            } else if (args_vm["compression-level"].as<int>() > MAX_COMPRESSION_LEVEL ) {
                compression_level = MAX_COMPRESSION_LEVEL;
            } else {
                compression_level = (uint8_t)args_vm["compression-level"].as<int>();
            }
        }
        
        if (args_vm.count("compression-algorithm") != 0) {
            std::string algstr = args_vm["compression-algorithm"].as<std::string>();
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
        
        if (args_vm.count("input-files") !=0 ) {
            fnames = args_vm["input-files"].as< std::vector<std::string> >();
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
    for (const auto &fname: fnames) {
        auto starttime = std::chrono::high_resolution_clock::now();
        std::string fname_hdf = fname;
        stringreplace(fname_hdf, std::filesystem::path(fname).extension(), ".h5");
        std::cout << fname << "\treading..." << std::flush;
        std::vector<long long> data;
        readTSVtags(fname, data);
        std::cout << "ok\twriting..." << std::flush;
        writeHDFtags(fname_hdf, data, compression_alg, compression_level);
        auto endtime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endtime - starttime).count();
        std::cout << "ok\t" << duration << "ms" << std::endl;
    }
    return 0;
}
