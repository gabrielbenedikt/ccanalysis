#include "tagsplit.h"

namespace bpo = boost::program_options;

int main(int argc, char **argv)
{ 
    uint16_t numparts = 0;
    uint32_t numtags = 0;
    /*
     * get files from prompt
     */
    std::vector<std::string> tagfiles;
    try {
        bpo::options_description args("Arguments");
        args.add_options()
            ("help,h", "prints this message")
            ("input-files", bpo::value<std::vector<std::string>>(), "Tag files to split")
            ("num-parts,p", bpo::value<uint16_t>(), "split file in <n> roughly equal parts")
            ("num-tags,t",  bpo::value<uint32_t>(), "split file in parts containing <n> tags each")
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
            std::cout << "Usage: ./tagsplit [-p|-t <arg>] input-files\n";
            std::cout << "Splits tag file in parts.\n";
            std::cout << args << "\n";
            return 0;
        }

        if (args_vm.count("input-files") != 0) {
            tagfiles = args_vm["input-files"].as< std::vector<std::string> >();
        } else {
            std::cout << "no input files specified." << std::endl;
            return 2;
        }
        
        if ((args_vm.count("num-parts") != 0) && (args_vm.count("num-tags") != 0)) {
            std::cout << "Error: incompatible split methods specified. Please pick exactly 1 method." << std::endl;
            return 2;
        }
        if ((args_vm.count("num-parts") == 0) && (args_vm.count("num-tags") == 0)) {
            std::cout << "Error: neither number of output files nor number of tags per file specified." << std::endl;
            return 2;
        }
        
        if (args_vm.count("num-parts") != 0) {
            numparts = args_vm["num-parts"].as< uint16_t >();
        }
        if (args_vm.count("num-tags") != 0) {
            numtags = args_vm["num-tags"].as< uint32_t >();
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
    
    for (const auto &fn: tagfiles) {
        std::cout << "splitting " << fn << std::endl;
        if ((std::filesystem::path(fn).extension() == ".h5")) {
            //TODO: implement
            std::cout << "HDF5 splitting not yet implemented" << std::endl;
            continue;
            split_hdf5(fn, numparts, numtags);
        } else if ((std::filesystem::path(fn).extension() == ".txt") || (std::filesystem::path(fn).extension() == ".tsv")) {
            split_tsv(fn, numparts, numtags);
        } else if ((std::filesystem::path(fn).extension() == ".tags") || (std::filesystem::path(fn).extension() == ".zst")) {
            //TODO: implement
            std::cout << "Capnproto splitting not yet implemented" << std::endl;
            continue;
            split_cap(fn, numparts, numtags);
        } else {
            std::cerr << "filetype not recognized." << std::endl;
            continue;
        }
    }
}

void split_tsv(const std::string &fn, const uint16_t in_numparts, const uint32_t in_numtags){
    if ((in_numparts==0)&&(in_numtags==1)) {
        std::cout << "Error: both number of parts and number of tags are 0." << std::endl;
        return;
    }
    
    uint64_t numlines = get_num_lines(fn);
    uint32_t numtags = in_numtags;
    std::string old_extension = std::filesystem::path(fn).extension();
    if (!(in_numparts==0)) {
        numtags = std::ceil(numlines/in_numparts);
    } 
    
    size_t numfiles = std::ceil(numlines/numtags);
    
    std::cout << "input file: " << fn << std::endl;
    std::cout << "#tags in input file / #tags per output file / #output files: " << numlines << " / " << numtags << " / " << numfiles << std::endl;
    
    if (numfiles == 1) {
        std::cout << "Number of tags in file is already smaller than requested number of tags per split file.\n";
        std::cout << "Doing nothing." << std::endl;
    } else {
        uint8_t counter_width = (uint8_t)std::floor(std::log10(numfiles)) +1;
        std::ifstream ifs(fn);
        std::stringstream ss;
        std::string s;
        ss << ifs.rdbuf();
        for (size_t i = 0; i<=numfiles; ++i) {
            std::cout << "writing output file " << i << " / " << numfiles << std::endl;
            std::string numstr = fixed_width_intstr(i, counter_width);
            std::string out_fn = std::filesystem::path(fn).replace_extension(numstr+old_extension);
            std::ofstream ofs(out_fn);
            uint64_t currentfilelines = 0;
            while(currentfilelines<=numtags) {
                std::getline(ss,s,'\n');
                if (s.empty()) {
                    break;
                }
                ofs << s << '\n';
                ++currentfilelines;
            }
            ofs.close();
        }
        ifs.close();
    }
}


void split_hdf5(const std::string &fn, const uint16_t numparts, const uint32_t numtags){
    if ((numparts==0)&&(numtags==1)) {
        std::cout << "Error: both number of parts and number of tags are 0." << std::endl;
        return;
    }
    
    if (!(numparts==0)) {
        ;;
    } 
    
    if (!(numtags==0)) {
        ;;
    }
    
    std::cout << "debug" << numparts << numtags << fn << std::endl;
}

void split_cap(const std::string &fn, const uint16_t numparts, const uint32_t numtags){
    if ((numparts==0)&&(numtags==1)) {
        std::cout << "Error: both number of parts and number of tags are 0." << std::endl;
        return;
    }
    
    if (!(numparts==0)) {
        ;;
    } 
    
    if (!(numtags==0)) {
        ;;
    }
    
    
    std::cout << "debug" << numparts << numtags << fn<< std::endl;
}
