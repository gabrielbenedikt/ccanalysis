#include "cap2tsv.h"

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
        bpo::options_description args_hidden("Hidden options");
        args_hidden.add_options()
            ("input-files", bpo::value<std::vector<std::string>>(), "HDF5 files to convert")
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
            std::cout << "Converts capnproto files to tab separated values (TSV) files.\n";
            std::cout << args << "\n";
            return 0;
        }

        if (args_vm.count("input-files") != 0) {
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
        std::string fn_tsv = fname;
        stringreplace(fn_tsv, std::filesystem::path(fname).extension(), ".txt");
        std::cout << fname << "\treading..." << std::flush;
        std::vector<long long> data;
        readcapnptags(fname, data);
        std::cout << "ok\twriting..." << std::flush;
        lltoTSV(fn_tsv, data);
        auto endtime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endtime - starttime).count();
        std::cout << "ok\t" << duration << "ms" << std::endl;
    }
    return 0;
}
