#include "histogram.h"

namespace bpo = boost::program_options;

int main(int argc, char **argv)
{   
    uint32_t ETA = 0;
    uint32_t SPENT_TIME = 0;
    
    Config cfg = read_config();

    omp_set_num_threads(cfg.NUM_THREADS);
    
    /*
     * get files from prompt
     */
    std::vector<std::string> newtagfiles;
    try {
        bpo::options_description args("Arguments");
        args.add_options()
            ("input-files", bpo::value<std::vector<std::string>>(), "TSV files to convert")
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
            std::cout << "Usage: ./histogram  input-files\n";
            std::cout << "creates histograms from tagfiles.\n";
            std::cout << args << "\n";
            return 0;
        }

        if (args_vm.count("input-files") != 0) {
            newtagfiles = args_vm["input-files"].as< std::vector<std::string> >();
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
    
    if (newtagfiles.empty()){
        newtagfiles = get_new_tagfiles(cfg);
    }

    //when multiple workers are iterating over the vector, this will prevent them from analyzing the same file
    std::random_device rd;
    std::mt19937 g(rd());
    std::ranges::shuffle ( newtagfiles, g);

    size_t nanalyzed = 0;
    for ( size_t i = 0; i<newtagfiles.size(); ++i ) {
        auto starttime = std::chrono::high_resolution_clock::now();
        std::string fn = newtagfiles[i];

        // check if there is already an analyzed datafile. could appear if several PCs are working on the same dataset
        std::string savefname = fn;
        std::string newext;
        if (cfg.TRUNCATE_S > 0) {
            std::string truncstr = std::string(3 - std::to_string(cfg.TRUNCATE_S).length(), '0') + std::to_string(cfg.TRUNCATE_S);
            newext = "_trunc"+ truncstr +"_ccs.pbdat";
        } else {
            newext = ".pbdat";
        }
        stringreplace(savefname,std::filesystem::path(savefname).extension(), newext);

        if (fileExists(savefname)) {
            std::cout << "file seems to have been already analyzed after list of datafiles was composed: \n";
            std::cout << savefname << std::endl;
            auto endtime = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(endtime - starttime).count();
            ++nanalyzed;
            SPENT_TIME += duration;
            continue;
        }

        std::cout << "file " << nanalyzed << " of " << newtagfiles.size() << ": filename: " << fn << std::endl;
        
        /*
        * read tag file
        */
        std::vector<long long> data;
        
        if ((std::filesystem::path(fn).extension() == ".h5")) {
            readHDF5tags(fn, data);
        } else if ((std::filesystem::path(fn).extension() == ".txt") || (std::filesystem::path(fn).extension() == ".tsv")) {
            readTSVtags(fn, data);
        } else if ((std::filesystem::path(fn).extension() == ".tags") || (std::filesystem::path(fn).extension() == ".zst")) {
            readcapnptags(fn, data);
        } else {
            std::cerr << "filetype not recognized." << std::endl;
            continue;
        }
        
        //get all channels in patterns
        std::vector<uint16_t> channels = unique(flatten(cfg.patterns));
        std::vector<std::vector<long long>> tags_per_channel = {};
        
        /*
        * separate into channels
        */
        separate_tags_per_channels(data, tags_per_channel, channels, cfg);

        //free memory of input data
        data = std::vector<long long>();
        
        /*
        * create histogram
        */
        std::vector<long long> offsets = arange<long long>(cfg.HIST_START,cfg.HIST_STOP,cfg.HIST_STEP);
        std::vector<histogram_onepattern> histvec = {};
        for (const auto &p : cfg.patterns) {
            histogram_onepattern hist = histogram(&tags_per_channel, &channels, &offsets, p, cfg.HIST_STEP);
            histvec.push_back(hist);
        }
        histograms allhistograms;
        histograms_to_struct(&histvec, &allhistograms, cfg);
        // save analysis to disk
        histstruct_protobuf_todisk(&allhistograms, savefname);
        
        // how much time did analysis take? give eta for all files
        auto endtime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endtime - starttime).count();
        ++nanalyzed;
        SPENT_TIME += duration; 
        ETA = (SPENT_TIME/nanalyzed)*(newtagfiles.size()-nanalyzed);
        
        std::cout << "analysis T / avg T / ETA: " << duration << "s / " << SPENT_TIME/(i+1) <<"s / " << ETA/(MIN_TO_SEC*omp_get_num_threads()) <<"m" << std::endl;
    }
    
    return 0;
}


/********************************************************************************
*** separate tag vectors
*/
void separate_tags_per_channels(const std::vector<long long> &tags, std::vector<std::vector<long long>>& tags_per_channel, const std::vector<uint16_t> &channels, const Config& cfg) {
    for (size_t i = 0; i<channels.size(); ++i) {
        tags_per_channel.emplace_back(std::vector<long long>());
    }
    
    uint64_t numtags = tags.size();
    uint64_t maxtag = numtags;
    if (cfg.TRUNCATE_S > 0) {
        long long firsttag = tags[1];
        long long lasttag = 0;
        for (uint64_t i=2; i<numtags; i+=2) {
            lasttag = tags[i+1];
            //some tags are zero. this breaks the comparison
            if (lasttag!=0) {
                if ((lasttag-firsttag)*cfg.CS*SEC_TO_NS>cfg.TRUNCATE_S) {
                    maxtag = i;
                    break;
                }
            }
        }
    }

    long long currtag = 0;
    for (uint64_t i=0; i<maxtag; i+=2) {
        currtag = tags[i+1];
        //some tags are zero. this breaks stuff. exclude them. effing tagger.
        if (currtag!=0) {
            auto it = std::find(channels.begin(), channels.end(), tags[i]);
            if (it != channels.end()) {
                size_t idx = it-channels.begin();
                tags_per_channel[idx].emplace_back(tags[i+1]);
            }
        }
    }

    /* make sure tags are ordered
     * 
     * may be more efficient if you do it directly on long long* tags, but doing it on the std::vectors is easier to write.
     * and it should happen quite rarely, so there's no big overhead to be expected
     * 
     * tags are, apart from an overflow, monotonically increasing.
     * the only thing we need to worry about is when tags change from positive to negative values
     * pos->pos is fine
     * neg->neg is fine (monotonically inceasing)
     * neg->pos is find (still monotonically inceasing)
     * pos->neg is bad.
     */
    if ((tags[1] > 0) && (tags[numtags-1] < 0)) {
        std::cout << "tagjump detected."  << std::endl;
        /* tagjump */
        // find discontinuity
        auto is_neg = [](long long x) { return ( x < 0 ); };
        const long long diff = 1LL<<42;                                   /////////////////////long long diff = *prev(it) - *it;
        for (std::vector<long long> v: tags_per_channel) {
            auto it = std::ranges::find_if(v, is_neg);
        
            // add the last positive element to all remaining elements in vector
            // by taking this difference, the first tag after the jump will be a 'double click'.
            // Better: find out range of tags and add this
            
            std::transform(it, v.end(), it, [](long long x) {return x+diff; } ); 
        }
        
    } // else same sign, thus no tag jump
}

/********************************************************************************
*** find coincidences - set inersect
*/
histogram_onepattern histogram(const std::vector<std::vector<long long>>* tags_per_channel, const std::vector<uint16_t>* channels, const std::vector<long long>* offsets, const std::vector<uint16_t>& pattern, const long long &wnd) {
    std::vector<long long> tags_trigger = tags_per_channel->at(std::find(channels->begin(), channels->end(), pattern[0])-channels->begin());
    //size_t trigger_idx = std::find(channels->begin(), channels->end(), pattern[0])-channels->begin();
    //auto &tags_trigger = tags_per_channel[trigger_idx];
    
    std::vector<long long> tags_idler = tags_per_channel->at(std::find(channels->begin(), channels->end(), pattern[1])-channels->begin());
    std::ranges::for_each(tags_trigger.begin(), tags_trigger.end(), [&wnd](long long &tt) { tt=roundto(tt, wnd); });

    size_t offset_len = offsets->size();
    
    histogram_onepattern histogram = {};
    histogram.pattern = pattern;
    histogram.meastime = tags_trigger.back()-tags_trigger.front();
    histogram.cc.resize(offset_len);
    histogram.cc_tags.resize(offset_len);
    histogram.offsets = *offsets;
    
    #pragma omp parallel for
    for(uint i=0; i<offset_len; ++i) {
        long long os = offsets->at(i);
        
        // H
        std::vector<long long> tmptags = tags_idler;
        std::ranges::for_each(tmptags.begin(), tmptags.end(), [&os, &wnd](long long &ht) { ht=roundto(ht-os,wnd); });
        std::vector<long long> cc_tags = {};
        cc_tags.reserve(4096);
        std::set_intersection(tags_trigger.begin(), tags_trigger.end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));
        
        histogram.cc[i] = cc_tags.size();
        histogram.cc_tags[i] = cc_tags;
    }
    
    return histogram;         
}

void histograms_to_struct(const std::vector<histogram_onepattern> *pts, histograms *hs, Config& cfg) {
    for (const auto &h: *pts) {
        hs->offsets.emplace_back(h.offsets);
        hs->cc.emplace_back(h.cc);
        hs->cc_tags.emplace_back(h.cc_tags);
        hs->pattern.emplace_back(h.pattern);
        hs->meastime.emplace_back(h.meastime);
    }
    hs->resolution = cfg.CS;
}

int histstruct_protobuf_todisk(const histograms* data, const std::string &fname) {
    std::cout << "writing to file " << fname << std::endl;
    histogramset::histograms hdat;
    for (size_t i = 0; i<data->meastime.size(); ++i) {
        hdat.add_meastime(data->meastime[i]);
        histogramset::histograms::i64arr* offsets = hdat.add_offsets();
        histogramset::histograms::i64arr* cc = hdat.add_cc();
        histogramset::histograms::u32arr* pattern = hdat.add_pattern();
        histogramset::histograms::repi64arr* cc_tags_vec = hdat.add_cc_tags();
        
        
        for (auto o: data->offsets[i]) {
            offsets->add_arr(o);
        }
        for (auto c: data->cc[i]) {
            cc->add_arr(c);
        }
        for (auto p: data->pattern[i]) {
            pattern->add_arr(p);
        }
        for (const auto &tarr: data->cc_tags[i]) {
            histogramset::histograms::i64arr* cc_tags = cc_tags_vec->add_arr();
            for (auto t: tarr) {
                    cc_tags->add_arr(t);
            }
        }
        /*
        for (size_t j = 0; j<data->cc_tags[i].size(); ++j) {
            histogramset::histograms::i64arr* cc_tags = cc_tags_vec->add_arr();
            for (size_t k = 0; k<data->cc_tags[i][j].size(); ++k) {
                    cc_tags->add_arr(data->cc_tags[i][j][k]);
            }
        }
        */
    }
    
    std::fstream ofs(fname, std::ios::out | std::ios::trunc | std::ios::binary);
    if (!hdat.SerializeToOstream(&ofs)) {
        std::cerr << "Failed to write." << std::endl;
        return -1;
    }
    ofs.close();
    
    return 0;
}


/********************************************************************************
*** find new tagfiles
*/
std::vector<std::string> get_new_tagfiles(const Config& cfg) {
    std::string path = std::filesystem::current_path();
    std::vector<std::string> alltagfiles = {};
    std::vector<std::string> analyzedfiles = {};
    std::vector<std::string> newtagfiles = {};
    
    std::string subtruncstr;
    if (cfg.TRUNCATE_S > 0) {
        std::string truncstr = std::string(3 - std::to_string(cfg.TRUNCATE_S).length(), '0') + std::to_string(cfg.TRUNCATE_S);
        subtruncstr = "_trunc"+ truncstr;
    } else {
        subtruncstr = "_trunc000";
    }
            
            
    for ( const auto & p : std::filesystem::directory_iterator(path) ) {
        if ( p.path().extension() == ".h5" ) {
            if (p.path().filename() != "tomos.h5") {
                alltagfiles.emplace_back(p.path());
            }
        } else if ( p.path().extension() == ".pbdat" ) {
            if (p.path().filename().string().find(subtruncstr) != std::string::npos) {
                analyzedfiles.emplace_back(p.path());
            } 
        }
    }
    
    newtagfiles = alltagfiles;
    if ( analyzedfiles.empty() ) {
        return alltagfiles;
    }
    
    auto it = newtagfiles.begin();
    while (it != newtagfiles.end() ) {
        bool found = false;
        std::string h5string = *it;
        stringreplace(h5string, ".h5","");
        for ( const auto &p : analyzedfiles ) {
            std::string tmpstring = p;
            stringreplace(tmpstring,subtruncstr+"_ccs.pbdat","");
            
            if (( tmpstring == h5string)) {
                it=newtagfiles.erase(it);
                found = true;
            }
        }
        if (!found) {
            ++it;
        }
    }
    
    return newtagfiles;
}

/********************************************************************************
*** read config file
*/
Config read_config() {
    Config cfg;
    std::ifstream cfgf("histogram.cfg");
    if (cfgf.is_open()) {
        std::string line;
        while(getline(cfgf, line)){
            line.erase(remove_if(line.begin(), line.end(), ::isspace),
                                 line.end());
            if(line[0] == '#' || line.empty()) {
                continue;
            }
            auto delimiterPos = line.find('=');
            auto name = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);
            
            if (name == "hist_start") {
                cfg.HIST_START = stoll(value);
            } else if (name == "hist_stop") {
                cfg.HIST_STOP = stoll(value);
            } else if (name == "hist_step") {
                cfg.HIST_STEP = stoll(value);
            } else if (name == "truncate") {
                cfg.TRUNCATE_S = stoul(value);
            } else if (name == "patterns") {
                cfg.patterns = parse_patterns(value);
            } else if (name == "num_threads") {
                cfg.NUM_THREADS = stoi(value);
            } else if (name == "tagger_resolution") {
                cfg.CS = stod(value);
            } else {
                std::cout << "unknown config name: " << name << std::endl;
            }
        }
        cfgf.close();
    } else {
        std::cerr << "Couldn't open config file for reading.\n";
    }

    std::cout << "---config---\n";
    std::cout << "hist_start    : " << cfg.HIST_START << "cs\n";
    std::cout << "hist_stop     : " << cfg.HIST_STOP  << "cs\n";
    std::cout << "hist_step     : " << cfg.HIST_STEP  << "cs\n";
    std::cout << "truncate at   : " << cfg.TRUNCATE_S << "s\n";
    std::cout << "tag resolution: " << cfg.CS         << "ns\n";
    std::cout << "# threads     : " << cfg.NUM_THREADS<< "\n";
    std::cout << "patterns:     : " << std::endl;
    for (const auto &e: cfg.patterns) {
        print_vector(e);
    }
    
    return cfg;
}
