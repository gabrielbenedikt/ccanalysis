#include "hdf2tsv.h"

namespace bpo = boost::program_options;

int main(int argc, char **argv) {
    
    std::vector<std::string> fnames = {};
    /*
     * parse arguments
     */
    try {
        bpo::options_description args("Arguments");
        args.add_options()
            ("help,h", "prints this message")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        bpo::options_description hidden("Hidden options");
        hidden.add_options()
            ("input-files", bpo::value<std::vector<std::string>>(), "HDF5 files to convert")
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
            std::cout << "Converts capnproto files to tab separated values (TSV) files.\n";
            std::cout << args << "\n";
            return 0;
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
        std::string fn_tsv = fn;
        stringreplace(fn_tsv, std::filesystem::path(fn).extension(), ".txt");
        std::cout << fn << "\treading..." << std::flush;
        long long cap_len = 0;
        std::vector<long long> data;
        readcapnptags(fn, data, cap_len);
        std::cout << "ok\twriting..." << std::flush;
        lltoTSV(fn_tsv, data, cap_len);
        auto endtime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endtime - starttime).count();
        std::cout << "ok\t" << duration << "ms" << std::endl;
    }
    return 0;
}
