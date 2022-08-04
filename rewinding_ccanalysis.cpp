#include "rewinding_ccanalysis.h"

int main(void)
{   
    read_config();

    omp_set_num_threads(cfg.NUM_THREADS);
    
    std::vector<std::string> newtagfiles = get_new_tagfiles();

    //when multiple workers are iterating over the vector, this will prevent them from doing the same thing
    std::random_device rd;
    std::mt19937 g(rd());
    std::ranges::shuffle ( newtagfiles, g);

    uint nanalyzed = 0;
    for ( uint i = 0; i<newtagfiles.size(); ++i ) {
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
        * read hdf file
        */
        long long data_len = 0;
        long long *data = readHDF5tags(fn, data_len);
        
        /*
        * separate into channels
        * already converted to ns
        */
        std::vector<double> cnt_tr = {};
        std::vector<double> cnt_h = {};
        std::vector<double> cnt_v = {};
        std::vector<double> cnt_fpga = {};
        separate_tags_per_channels(data, data_len, cnt_tr, cnt_h, cnt_v, cnt_fpga);
        
        std::ranges::for_each(cnt_tr.begin(), cnt_tr.end(), [](double &tt) { tt=round(tt*10)/10; });

        /*
        * find coincidences
        */
        double meastime = (cnt_tr.back()-cnt_tr.front())*pow(10,-9);
        std::vector<double> offsets = arange<double>(cfg.HIST_START,cfg.HIST_STOP,cfg.HIST_STEP);
        std::vector<cc_point> hist_points;
        if (cfg.FPGA_USE > 0) {
            std::cout << "with fpga" << std::endl;
            hist_points = find_coincidences_fpga(&cnt_tr, &cnt_h, &cnt_v, &offsets, &cnt_fpga);
        } else {
            std::cout << "wo fpga" << std::endl;
            hist_points = find_coincidences(&cnt_tr, &cnt_h, &cnt_v, &offsets);
        }

        //convert point vector to single struct
        ccstruct hist_combined;
        cc_point_to_ccstruct(&hist_points, &hist_combined);
        hist_combined.meastime = meastime;
        hist_combined.num_trigger_tags = cnt_tr.size();
        hist_combined.num_fpga_tags = cnt_fpga.size();

        // save analysis to disk
        ccstruct_protobuf_todisk(&hist_combined, savefname);

        // how much time did analysis take? give eta for all files
        auto endtime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endtime - starttime).count();
        ++nanalyzed;
        SPENT_TIME += duration; 
        ETA = (SPENT_TIME/nanalyzed)*(newtagfiles.size()-nanalyzed);
        
        std::cout << "number of trigger / H / V tags: " << cnt_tr.size() << " / " << cnt_h.size() << " / " << cnt_v.size() << "\n";
        std::cout << "measurement time: " << meastime << "\n";
        std::cout << "analysis T / avg T / ETA: " << duration << "s / " << SPENT_TIME/i <<"s / " << ETA/(60*omp_get_num_threads()) <<"m" << std::endl;
        
        delete[] data;
    }
    
    std::cout << "done" << std::endl;
    
    return 0;
}


/********************************************************************************
*** separate tag vectors
*/
void separate_tags_per_channels(const long long* tags, const long long numtags, std::vector<double>& cnt_tr, std::vector<double>& cnt_h, std::vector<double>& cnt_v, std::vector<double>& cnt_fpga)
{   
    cnt_tr.reserve(numtags);
    cnt_fpga.reserve(numtags);
    cnt_h.reserve(numtags);
    cnt_v.reserve(numtags);
    
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
            if (tags[i]==cfg.CHN_TR) {
                cnt_tr.emplace_back(tags[i+1]*cfg.CS);
            } else if (tags[i]==cfg.CHN_H) {
                cnt_h.emplace_back(tags[i+1]*cfg.CS);
            } else if (tags[i]==cfg.CHN_V) {
                cnt_v.emplace_back(tags[i+1]*cfg.CS);
            } else if (tags[i]==cfg.CHN_FPGA) {
                cnt_fpga.emplace_back(tags[i+1]*cfg.CS);
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
    if ((cnt_tr.front() > 0) && (cnt_tr.back() < 0)) {
        std::cout << "tagjump detected."  << std::endl;
        /* tagjump */
        // find discontinuity
        auto is_neg = [](double x) { return ( x < 0 ); };
        auto it = std::ranges::find_if(cnt_tr, is_neg);
        
        // add the last positive element to all remaining elements in vector
        // by taking this difference, the first tag after the jump will be a 'double click'.
        // Better: find out range of tags and add this
        
        double diff = pow(2,42);                                    /////////////////////double diff = *prev(it) - *it;
        std::transform(it, cnt_tr.end(), it, [&diff](double x) {return x+diff; } ); 
        
        //now, do the same for cnt_h and cnt_v
        if ((cnt_h.front() > 0) && (cnt_h.back() < 0)) {
            auto it = std::ranges::find_if(cnt_h, is_neg);
            //diff = *prev(it) - *it;
            std::transform(it, cnt_h.end(), it, [&diff](double x) {return x+diff; } ); 
        }
        if ((cnt_v.front() > 0) && (cnt_v.back() < 0)) {
            auto it = std::ranges::find_if(cnt_v, is_neg);
            //diff = *prev(it) - *it;
            std::transform(it, cnt_v.end(), it, [&diff](double x) {return x+diff; } ); 
        }
        if ((cnt_fpga.front() > 0) && (cnt_fpga.back() < 0)) {
            auto it = std::ranges::find_if(cnt_fpga, is_neg);
            //diff = *prev(it) - *it;
            std::transform(it, cnt_fpga.end(), it, [&diff](double x) {return x+diff; } ); 
        }
    } // else same sign, thus no tag jump
}



/********************************************************************************
*** find coincidences - set inersect
*/
std::vector<cc_point> find_coincidences(const std::vector<double>* ttags, const std::vector<double>* htags, const std::vector<double>* vtags, const std::vector<double>* offsets) { 
    const uint offset_len = offsets->size();
    std::vector<cc_point> result(offset_len);
    
    #pragma omp parallel for
    for(uint i=0; i<offset_len; ++i) {
        double os = offsets->at(i);
        
        cc_point ccdata;
        
        // H
        std::vector<double> tmptags = *htags;
        std::ranges::for_each(tmptags.begin(), tmptags.end(), [&os](double &ht) { ht= std::round((ht-os)*10)/10; });
        std::vector<double> cc_tags = {};
        cc_tags.reserve(4096);
        std::set_intersection(ttags->begin(), ttags->end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));

        ccdata.cc_h = cc_tags.size();
        ccdata.cc_h_tags = cc_tags;
        
        // V
        tmptags = *vtags;
        std::ranges::for_each(tmptags.begin(), tmptags.end(), [&os](double &vt) { vt = std::round((vt-os)*10)/10; });
        
        cc_tags = {};
        std::set_intersection(ttags->begin(), ttags->end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));



        ccdata.cc_v = cc_tags.size();
        ccdata.cc_v_tags = cc_tags;
        ccdata.offset = os;
        result[i] = ccdata;
    }
    return result;         
}

/********************************************************************************
*** find coincidences with trigger - set inersect - keep only coincidences where fpga signal is present too
*/
std::vector<cc_point> find_coincidences_fpga(const std::vector<double>* ttags, const std::vector<double>* htags, const std::vector<double>* vtags, const std::vector<double>* offsets, const std::vector<double>* ftags) {
    //round fpga tags to 10 ns
    std::vector<double> fpga_tags_rounded = *ftags;
    std::ranges::for_each(fpga_tags_rounded.begin(), fpga_tags_rounded.end(), [](double &ft) { ft= std::round(ft); });

    const uint offset_len = offsets->size();
    std::vector<cc_point> result(offset_len);

    #pragma omp parallel for
    for(uint i=0; i<offset_len; ++i) {
        double os = offsets->at(i);

        cc_point ccdata;

        // H
        std::vector<double> tmptags = *htags;
        std::ranges::for_each(tmptags.begin(), tmptags.end(), [&os](double &ht) { ht= std::round((ht-os)*10)/10; });
        std::vector<double> cc_tags = {};
        cc_tags.reserve(4096);
        std::set_intersection(ttags->begin(), ttags->end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));

        //check coincidences between fpga and cc_tags
        std::vector<double> ccc_tags = {};
        for (double cctag : cc_tags) {
            auto lbh = std::lower_bound(fpga_tags_rounded.begin(), fpga_tags_rounded.end(), cctag+cfg.FPGA_DELAY_MIN);
            auto ubh = std::upper_bound(fpga_tags_rounded.begin(), fpga_tags_rounded.end(), cctag+cfg.FPGA_DELAY_MAX);
            if (lbh == ubh || lbh == fpga_tags_rounded.end() || ubh == fpga_tags_rounded.end()) {
            } else {
                ccc_tags.emplace_back(cctag);
            }
        }
        ccdata.cc_h = ccc_tags.size();
        ccdata.cc_h_tags = ccc_tags;


        // V
        tmptags = *vtags;
        std::ranges::for_each(tmptags.begin(), tmptags.end(), [&os](double &vt) { vt = std::round((vt-os)*10)/10; });

        cc_tags = {};
        std::set_intersection(ttags->begin(), ttags->end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));

        //check coincidences between fpga and cc_tags
        ccc_tags = {};
        ccc_tags.reserve(4096);
        for (double cctag : cc_tags) {
            auto lbv = std::lower_bound(fpga_tags_rounded.begin(), fpga_tags_rounded.end(), cctag+cfg.FPGA_DELAY_MIN);
            auto ubv = std::upper_bound(fpga_tags_rounded.begin(), fpga_tags_rounded.end(), cctag+cfg.FPGA_DELAY_MAX);
            if (lbv == ubv || lbv == fpga_tags_rounded.end() || ubv == fpga_tags_rounded.end()) {
            } else {
                ccc_tags.emplace_back(cctag);
            }
        }
        ccdata.cc_v = ccc_tags.size();
        ccdata.cc_v_tags = ccc_tags;

        ccdata.offset = os;
        result[i] = ccdata;
    }
    return result;
}

/********************************************************************************
*** separate tag vectors
*/
void cc_point_to_ccstruct(const std::vector<cc_point> *pts, ccstruct *ccs) {
    uint numpts = pts->size();
    std::vector<long> cc_h = {};
    cc_h.reserve(numpts);
    std::vector<long> cc_v = {};
    cc_v.reserve(numpts);
    std::vector<std::vector<double>> cc_h_tags = {{}};
    cc_h_tags.reserve(numpts);
    std::vector<std::vector<double>> cc_v_tags = {{}};
    cc_v_tags.reserve(numpts);
    std::vector<double> offsets = {};
    offsets.reserve(numpts);
    for (auto pnt: *pts) {
        cc_h.emplace_back(pnt.cc_h);
        cc_v.emplace_back(pnt.cc_v);
        cc_h_tags.emplace_back(pnt.cc_h_tags);
        cc_v_tags.emplace_back(pnt.cc_v_tags);
        offsets.emplace_back(pnt.offset);
    }
    
    ccs->offsets = offsets;
    ccs->cc_h = cc_h;
    ccs->cc_v = cc_v;
    ccs->cc_h_tags = cc_h_tags;
    ccs->cc_v_tags = cc_v_tags;
}

/********************************************************************************
*** print contents of a ccstruct
*/
void print_ccstruct(const ccstruct* dat) {
    std::cout << "offsets" << std::endl;
    for (uint i=0; i<dat->offsets.size(); ++i) {
        std::cout << dat->offsets[i] << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "cc_h" << std::endl;
    for (uint i=0; i<dat->cc_h.size(); ++i) {
        std::cout << dat->cc_h[i] << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "cc_v" << std::endl;
    for (uint i=0; i<dat->cc_v.size(); ++i) {
        std::cout << dat->cc_v[i] << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "cc_h_tags" << std::endl;
    for (uint i=0; i<dat->cc_h_tags.size(); ++i) {
        for (uint j=0; j<dat->cc_h_tags[i].size(); ++j) {
            std::cout << dat->cc_h_tags[i][j] << ", ";
        }
        std::cout << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "cc_v_tags" << std::endl;
    for (uint i=0; i<dat->cc_v_tags.size(); ++i) {
        for (uint j=0; j<dat->cc_v_tags[i].size(); ++j) {
            std::cout << dat->cc_v_tags[i][j] << ", ";
        }
         std::cout << ", ";
    }
    
    std::cout << "meastime: " << dat->meastime << "\n";
    std::cout << "number of trigger tags: " << dat->num_trigger_tags << "\n";
    std::cout << "number of fpga tags: " << dat->num_fpga_tags << "\n";
    std::cout << std::endl;
}

/********************************************************************************
*** save a ccstruct to disk using protobuf
*/
int ccstruct_protobuf_todisk(const ccstruct* data, const std::string fname) {
    std::cout << "writing to file " << fname << std::endl;
    ccset::ccset_data pbdat;
    
    for (uint i=0; i< data->offsets.size(); ++i) {
        pbdat.add_offsets(data->offsets[i]);
        pbdat.add_cc_h(data->cc_h[i]);
        pbdat.add_cc_v(data->cc_v[i]);
        ccset::ccset_data::darr* arr2dh = pbdat.add_cc_h_tags();
        ccset::ccset_data::darr* arr2dv = pbdat.add_cc_v_tags();

        for (uint j=0; j< data->cc_h_tags[i].size(); ++j) {
            arr2dh->add_arr(data->cc_h_tags[i][j]);
        }
        for (uint j=0; j< data->cc_v_tags[i].size(); ++j) {
            arr2dv->add_arr(data->cc_v_tags[i][j]);
        }

    }
    pbdat.set_meastime(data->meastime);
    pbdat.set_num_trigger_tags(data->num_trigger_tags);
    pbdat.set_num_fpga_tags(data->num_fpga_tags);
    
    /*
     * write to disk
     */
    {
        std::fstream ofs(fname, std::ios::out | std::ios::trunc | std::ios::binary);
        if (!pbdat.SerializeToOstream(&ofs)) {
            std::cerr << "Failed to write." << std::endl;
            return -1;
        }
        ofs.close();
    }
    
    return 0;
}


/********************************************************************************
*** load a ccstruct from disk using protobuf
*/
int ccstruct_protobuf_fromdisk(ccstruct* data, const std::string fname)
{
    //read from disk
    ccset::ccset_data read_pbdat;
    {
        std::fstream ifs(fname, std::ios::in | std::ios::binary);
        if (!read_pbdat.ParseFromIstream(&ifs)) {
            std::cerr << "Failed to parse." << std::endl;
            return -1;
        }
        ifs.close();
    }
    
    // protobuf data to struct
    data->offsets = {};
    data->cc_h = {};
    data->cc_v = {};
    data->cc_h_tags = {  };
    data->cc_v_tags = {  };
    for (int i=0; i<read_pbdat.offsets_size(); ++i) {
        data->offsets.emplace_back(read_pbdat.offsets(i));
    }
    for (int i=0; i<read_pbdat.cc_h_size(); ++i) {
        data->cc_h.emplace_back(read_pbdat.cc_h(i));
    }
    for (int i=0; i<read_pbdat.cc_v_size(); ++i) {
        data->cc_v.emplace_back(read_pbdat.cc_v(i));
    }
    for (int i=0; i<read_pbdat.cc_h_tags_size(); ++i) {
        std::vector<double> tmparr = { };
        for (int j=0; j<read_pbdat.cc_h_tags(i).arr_size(); ++j) {
            tmparr.emplace_back(read_pbdat.cc_h_tags(i).arr(j));
        }
        data->cc_h_tags.emplace_back(tmparr);
    }
    for (int i=0; i<read_pbdat.cc_v_tags_size(); ++i) {
        std::vector<double> tmparr = { };
        for (int j=0; j<read_pbdat.cc_v_tags(i).arr_size(); ++j) {
            tmparr.emplace_back(read_pbdat.cc_v_tags(i).arr(j));
        }
        data->cc_v_tags.emplace_back(tmparr);
    }
    data->meastime=read_pbdat.meastime();
    data->num_trigger_tags=read_pbdat.num_trigger_tags();
    data->num_fpga_tags=read_pbdat.num_fpga_tags();
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
    std::ifstream cfgf("cc.cfg");
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
            } else if (name == "chn_h") {
                cfg.CHN_H = stoi(value);
            } else if (name == "chn_v") {
                cfg.CHN_V = stoi(value);
            } else if (name == "chn_tr") {
                cfg.CHN_TR = stoi(value);
            } else if (name == "chn_fpga") {
                cfg.CHN_FPGA = stoi(value);
            } else if (name == "use_fpga") {
                cfg.FPGA_USE = stoi(value);
            } else if (name == "fpga_offset") {
                cfg.FPGA_OFFSET = stod(value);
            } else if (name == "fpga_delay_min") {
                cfg.FPGA_DELAY_MIN = stod(value);
            } else if (name == "fpga_delay_max") {
                cfg.FPGA_DELAY_MAX = stod(value);
            } else if (name == "fpga_hist_step") {
                cfg.FPGA_HIST_STEP = stod(value);
            } else if (name == "truncate") {
                cfg.TRUNCATE_S = stoi(value);
            } else if (name == "patterns") {
                cfg.patterns = parse_patterns(value);
            }else if ("num_threads") {
                cfg.NUM_THREADS = stoi(value);
            }else if ("tagger_resolution") {
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
    std::cout << "thigger chn   : " << cfg.CHN_TR     << "\n";
    std::cout << "fpga chn      : " << cfg.CHN_FPGA   << "\n";
    std::cout << "H chn         : " << cfg.CHN_H      << "\n";
    std::cout << "V chn         : " << cfg.CHN_V      << "\n";
    std::cout << "truncate at   : " << cfg.TRUNCATE_S << "s\n";
    std::cout << "use fpga chan : " << cfg.FPGA_USE   << "\n";
    std::cout << "fpga delay min: " << cfg.FPGA_DELAY_MIN << "ns\n";
    std::cout << "fpga delay max: " << cfg.FPGA_DELAY_MAX << "ns\n";
    std::cout << "patterns:     : " << std::endl;
    for (auto e: cfg.patterns) {
        print_vector(e);
    }

}



