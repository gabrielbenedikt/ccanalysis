#include "histogram.h"

int main(void)
{   
    read_config();

    omp_set_num_threads(cfg.NUM_THREADS);
    
    std::vector<std::string> newtagfiles = get_new_tagfiles();

    //when multiple workers are iterating over the vector, this will prevent them from doing the same thing
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
            newext = "_trunc000_ccs.pbdat";
        }
        stringreplace(savefname,".h5", newext);

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
        long long data_len = 0;
        long long *data;
        
        if ((std::filesystem::path(fn).extension() == ".h5")) {
            data = readHDF5tags(fn, data_len);
        } else if ((std::filesystem::path(fn).extension() == ".txt") || (std::filesystem::path(fn).extension() == ".tsv")) {
            std::vector<long long> r;
            readTSVtags(fn, r, data_len);
            data = &r[0];
        } else if ((std::filesystem::path(fn).extension() == ".tags") || (std::filesystem::path(fn).extension() == ".zst")) {
            std::vector<long long> r;
            r = readcapnptags(fn, data_len);
            data = &r[0];
        } else {
            std::cerr << "filetype not recognized." << std::endl;
            continue;
        }
        
        
        //get all channels in patterns
        std::vector<uint16_t> channels = unique(flatten(cfg.patterns));
        std::vector<std::vector<double>> tags_per_channel = {};
        /*
        * separate into channels
        */
        
        separate_tags_per_channels(data, data_len, tags_per_channel, channels);
        
        /*
        * create histogram
        */
        std::vector<double> offsets = arange<double>(cfg.HIST_START,cfg.HIST_STOP,cfg.HIST_STEP);
        std::vector<histogram_onepattern> histvec = {};
        for (auto p : cfg.patterns) {
            histogram_onepattern hist = histogram(&tags_per_channel, &channels, &offsets, p);
            histvec.push_back(hist);
        }
        histograms allhistograms;
        histograms_to_struct(&histvec, &allhistograms);
        // save analysis to disk
        histstruct_protobuf_todisk(&allhistograms, savefname);

        
        // how much time did analysis take? give eta for all files
        auto endtime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endtime - starttime).count();
        ++nanalyzed;
        SPENT_TIME += duration; 
        ETA = (SPENT_TIME/nanalyzed)*(newtagfiles.size()-nanalyzed);
        
         std::cout << "analysis T / avg T / ETA: " << duration << "s / " << SPENT_TIME/i <<"s / " << ETA/(60*omp_get_num_threads()) <<"m" << std::endl;
        
        delete[] data;
    }
    
    std::cout << "done" << std::endl;
    
    return 0;
}


/********************************************************************************
*** separate tag vectors
*/
void separate_tags_per_channels(const long long* tags, const long long numtags, std::vector<std::vector<double>>& tags_per_channel, const std::vector<uint16_t> channels) 
{   
    //TODO
    for (size_t i = 0; i<channels.size(); ++i) {
        tags_per_channel.emplace_back(std::vector<double>());
    }
    
    long long maxtag = numtags;
    if (cfg.TRUNCATE_S > 0) {
        double firsttag = tags[1]*cfg.CS;
        double lasttag = 0;
        for (long long i=2; i<numtags; i+=2) {
            lasttag = tags[i+1]*cfg.CS;
            //some tags are zero. this breaks the comparison
            if (lasttag!=0) {
                if ((lasttag-firsttag)>cfg.TRUNCATE_S*pow(10,9)) {
                    maxtag = i;
                    break;
                }
            }
        }
    }

    double currtag;
    for (long long i=0; i<maxtag; i+=2) {
        currtag = tags[i+1]*cfg.CS;
        //some tags are zero. this breaks stuff. exclude them. effing tagger.
        if (currtag!=0) {
            auto it = std::find(channels.begin(), channels.end(), tags[i]);
            if (it != channels.end()) {
                size_t idx = it-channels.begin();
                tags_per_channel[idx].emplace_back(tags[i+1]*cfg.CS);
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
        auto is_neg = [](double x) { return ( x < 0 ); };
        double diff = pow(2,42);                                    /////////////////////double diff = *prev(it) - *it;
        for (std::vector<double> v: tags_per_channel) {
            auto it = std::ranges::find_if(v, is_neg);
        
            // add the last positive element to all remaining elements in vector
            // by taking this difference, the first tag after the jump will be a 'double click'.
            // Better: find out range of tags and add this
            
            std::transform(it, v.end(), it, [&diff](double x) {return x+diff; } ); 
        }
        
    } // else same sign, thus no tag jump
}



/********************************************************************************
*** find coincidences - set inersect
*/
histogram_onepattern histogram(const std::vector<std::vector<double>>* tags_per_channel, const std::vector<uint16_t>* channels, const std::vector<double>* offsets, const std::vector<uint16_t> pattern) {
    std::vector<double> tags_trigger = tags_per_channel->at(std::find(channels->begin(), channels->end(), pattern[0])-channels->begin());
    std::vector<double> tags_idler = tags_per_channel->at(std::find(channels->begin(), channels->end(), pattern[1])-channels->begin());
    std::ranges::for_each(tags_trigger.begin(), tags_trigger.end(), [](double &tt) { tt=round(tt*10)/10; });
    
    histogram_onepattern histogram;
    histogram.pattern = pattern;
    histogram.meastime = tags_trigger.back()-tags_trigger.front();
    histogram.offsets = *offsets;
    
    #pragma omp parallel for
    for(uint i=0; i<offsets->size(); ++i) {
        double os = offsets->at(i);
        
        // H
        std::vector<double> tmptags = tags_idler;
        std::ranges::for_each(tmptags.begin(), tmptags.end(), [&os](double &ht) { ht= std::round((ht-os)*10)/10; });
        std::vector<double> cc_tags = {};
        cc_tags.reserve(4096);
        std::set_intersection(tags_trigger.begin(), tags_trigger.end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));

        
        histogram.cc.emplace_back(cc_tags.size());
        histogram.cc_tags.emplace_back(cc_tags);
        
    }
    return histogram;         
}

void histograms_to_struct(const std::vector<histogram_onepattern> *pts, histograms *hs) {
    for (auto h: *pts) {
        hs->offsets.emplace_back(h.offsets);
        hs->cc.emplace_back(h.cc);
        hs->cc_tags.emplace_back(h.cc_tags);
        hs->pattern.emplace_back(h.pattern);
        hs->meastime.emplace_back(h.meastime);
    }
}

int histstruct_protobuf_todisk(const histograms* data, const std::string fname) {
    std::cout << "writing to file " << fname << std::endl;
    histogramset::histograms hdat;
    for (size_t i = 0; i<data->meastime.size(); ++i) {
        hdat.add_meastime(data->meastime[i]);
        histogramset::histograms::darr* offsets = hdat.add_offsets();
        histogramset::histograms::i64arr* cc = hdat.add_cc();
        histogramset::histograms::u32arr* pattern = hdat.add_pattern();
        histogramset::histograms::repdarr* cc_tags_vec = hdat.add_cc_tags();
        
        for (size_t j = 0; j<data->offsets[i].size(); ++j) {
            offsets->add_arr(data->offsets[i][j]);
        }
        for (size_t j = 0; j<data->cc[i].size(); ++j) {
            cc->add_arr(data->cc[i][j]);
        }
        for (size_t j = 0; j<data->pattern[i].size(); ++j) {
            pattern->add_arr(data->pattern[i][j]);
        }
        for (size_t j = 0; j<data->cc_tags[i].size(); ++j) {
            histogramset::histograms::darr* cc_tags = cc_tags_vec->add_arr();
            for (size_t k = 0; data->cc_tags[i][j].size(); ++k) {
                cc_tags->add_arr(data->cc_tags[i][j][k]);
            }
        }
    }
    
    
    {
        std::fstream ofs(fname, std::ios::out | std::ios::trunc | std::ios::binary);
        if (!hdat.SerializeToOstream(&ofs)) {
            std::cerr << "Failed to write." << std::endl;
            return -1;
        }
        ofs.close();
    }
    
    return 0;
}


/********************************************************************************
*** find new tagfiles
*/
std::vector<std::string> get_new_tagfiles() {
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
            
            
    for ( auto & p : std::filesystem::directory_iterator(path) ) {
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
    if ( analyzedfiles.size() == 0 ) {
        return alltagfiles;
    } else {
        auto it = newtagfiles.begin();
        while (it != newtagfiles.end() ) {
            bool found = false;
            std::string h5string = *it;
            stringreplace(h5string, ".h5","");
            for ( auto p : analyzedfiles ) {
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
        
//         newtagfiles.erase(remove(newtagfiles.begin(), newtagfiles.end(), "tomos.h5"), newtagfiles.end());
        
        return newtagfiles;
    }
}

/********************************************************************************
*** read config file
*/
void read_config() {
    std::ifstream cfgf("histogram.cfg");
    if (cfgf.is_open()) {
        std::string line;
        while(getline(cfgf, line)){
            line.erase(remove_if(line.begin(), line.end(), ::isspace),
                                 line.end());
            if(line[0] == '#' || line.empty())
                continue;
            auto delimiterPos = line.find("=");
            auto name = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);
            
            if (name == "hist_start") {
                cfg.HIST_START = stod(value);
            } else if (name == "hist_stop") {
                cfg.HIST_STOP = stod(value);
            } else if (name == "hist_step") {
                cfg.HIST_STEP = stod(value);
            } else if (name == "truncate") {
                cfg.TRUNCATE_S = stoi(value);
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
    std::cout << "hist_start    : " << cfg.HIST_START << "ns\n";
    std::cout << "hist_stop     : " << cfg.HIST_STOP  << "ns\n";
    std::cout << "hist_step     : " << cfg.HIST_STEP  << "ns\n";
    std::cout << "truncate at   : " << cfg.TRUNCATE_S << "s\n";
    std::cout << "patterns:     : " << std::endl;
    std::cout << "tag resolution: " << cfg.CS           << "\n";
    std::cout << "# threads     : " << cfg.NUM_THREADS  << "\n";
    for (auto e: cfg.patterns) {
        print_vector(e);
    }
}
