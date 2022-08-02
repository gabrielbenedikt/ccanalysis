#include "rewinding_ccanalysis.h"

using namespace std;
//using namespace H5;

int main(void)
{   
    read_config();

    omp_set_num_threads(cfg.NUM_THREADS);
    
    vector<string> newtagfiles = get_new_tagfiles();

    //when multiple workers are iterating over the vector, this will prevent them from doing the same thing
    std::random_device rd;
    std::mt19937 g(rd());
    std::ranges::shuffle ( newtagfiles, g);

    uint nanalyzed = 0;
    for ( uint i = 0; i<newtagfiles.size(); ++i ) {
        auto starttime = std::chrono::high_resolution_clock::now();
        string fn = newtagfiles[i];

        // check if there is already an analyzed datafile. could appear if several PCs are working on the same dataset
        string savefname = fn;
        string newext;
        if (cfg.TRUNCATE_S > 0) {
            string truncstr = string(3 - to_string(cfg.TRUNCATE_S).length(), '0') + to_string(cfg.TRUNCATE_S);
            newext = "_trunc"+ truncstr +"_ccs.pbdat";
        } else {
            newext = "_trunc000_ccs.pbdat";
        }
        stringreplace(savefname,".h5", newext);

        if (fileExists(savefname)) {
            cout << "file seems to have been already analyzed after list of datafiles was composed: " << endl;
            cout << savefname << endl;
            auto endtime = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(endtime - starttime).count();
            ++nanalyzed;
            SPENT_TIME += duration;
            continue;
        }

        cout << "file " << nanalyzed << " of " << newtagfiles.size() << ": filename: " << fn << endl;
        
        /*
        * read hdf file
        */
        long long data_len = 0;
        long long *data = readHDF5tags(fn, data_len);
        
        /*
        * separate into channels
        * already converted to ns
        */
        vector<double> cnt_tr = {};
        vector<double> cnt_h = {};
        vector<double> cnt_v = {};
        vector<double> cnt_fpga = {};
        separate_tags_per_channels(data, data_len, cnt_tr, cnt_h, cnt_v, cnt_fpga);
        
        std::ranges::for_each(cnt_tr.begin(), cnt_tr.end(), [](double &tt) { tt=round(tt*10)/10; });

        /*
        * find coincidences
        */
        double meastime = (cnt_tr.back()-cnt_tr.front())*pow(10,-9);
        vector<double> offsets = arange<double>(cfg.HIST_START,cfg.HIST_STOP,cfg.HIST_STEP);
        
        vector<cc_point> hist_points;
        if (cfg.FPGA_USE > 0) {
            cout << "with fpga" << endl;
            hist_points = find_coincidences_fpga(&cnt_tr, &cnt_h, &cnt_v, &offsets, &cnt_fpga);
        } else {
            cout << "wo fpga" << endl;
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
        
        cout << "number of trigger / H / V tags: " << cnt_tr.size() << " / " << cnt_h.size() << " / " << cnt_v.size() << endl;
        cout << "measurement time: " << meastime << endl;
        cout << "analysis T / avg T / ETA: " << duration << "s / " << SPENT_TIME/i <<"s / " << ETA/(60*omp_get_num_threads()) <<"m" << endl;
        
        delete[] data;
    }
    
    cout << "done" << endl;
    
    return 0;
}

/********************************************************************************
*** read hdf5 tag file
*/
long long* readHDF5tags(const string fn, long long& out_data_len)
{
    string datasetPath = "/tags/block0_values";
    
    // Open HDF5 file handle, read only
    H5::H5File file(fn.c_str(), H5F_ACC_RDONLY);

    // access the required dataset by path name
    H5::DataSet dataset = file.openDataSet(datasetPath.c_str());

    // get the dataspace
    H5::DataSpace dataspace = dataset.getSpace();
  
    /*
    * Get the dimension size of each dimension in the dataspace and
    * display them.
    */
    hsize_t dims_out[2];
    /*int dimensions = */ dataspace.getSimpleExtentDims( dims_out, NULL);

    int NX=dims_out[0];
    int NY=dims_out[1];

    assert(NY==2);
    /*
    * Output buffer initialization.
    */
    long long* data_out = new long long[NX*NY];
    
    out_data_len = NX*NY;
    /*
    * Read data from hyperslab in the file into the hyperslab in
    * memory and display the data.
    */
    dataset.read( data_out, H5::PredType::NATIVE_LLONG);
    
    // close the HDF5 file
    file.close();
    
    return data_out;
}

std::vector<long long> readTSVtags(const std::string fn, long long& out_data_len)
{
    std::ifstream tsvfile(fn);
    std::string line;
    std::vector<long long> result;
    long long value;
    
    if (!tsvfile.is_open()) {
        throw std::runtime_error("Could not open file " + fn);
    }
    
    while (std::getline(tsvfile, line)) {
        std::stringstream ss(line);
        while (ss >> value) {
            result.push_back(value);
            if(ss.peek() == '\t') {
                ss.ignore();
            }
        }
    }
    
    out_data_len = result.size();
    
    return result;
}

/********************************************************************************
*** separate tag vectors
*/
void separate_tags_per_channels(const long long* tags, const long long numtags, vector<double>& cnt_tr, vector<double>& cnt_h, vector<double>& cnt_v, vector<double>& cnt_fpga)
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
        cout << "tagjump detected."  << endl;
        /* tagjump */
        // find discontinuity
        auto is_neg = [](double x) { return ( x < 0 ); };
        auto it = ranges::find_if(cnt_tr, is_neg);
        
        // add the last positive element to all remaining elements in vector
        // by taking this difference, the first tag after the jump will be a 'double click'.
        // Better: find out range of tags and add this
        
        double diff = pow(2,42);                                    /////////////////////double diff = *prev(it) - *it;
        transform(it, cnt_tr.end(), it, [&diff](double x) {return x+diff; } ); 
        
        //now, do the same for cnt_h and cnt_v
        if ((cnt_h.front() > 0) && (cnt_h.back() < 0)) {
            auto it = ranges::find_if(cnt_h, is_neg);
            //diff = *prev(it) - *it;
            transform(it, cnt_h.end(), it, [&diff](double x) {return x+diff; } ); 
        }
        if ((cnt_v.front() > 0) && (cnt_v.back() < 0)) {
            auto it = ranges::find_if(cnt_v, is_neg);
            //diff = *prev(it) - *it;
            transform(it, cnt_v.end(), it, [&diff](double x) {return x+diff; } ); 
        }
        if ((cnt_fpga.front() > 0) && (cnt_fpga.back() < 0)) {
            auto it = ranges::find_if(cnt_fpga, is_neg);
            //diff = *prev(it) - *it;
            transform(it, cnt_fpga.end(), it, [&diff](double x) {return x+diff; } ); 
        }
    } // else same sign, thus no tag jump
}



/********************************************************************************
*** find coincidences - set inersect
*/
vector<cc_point> find_coincidences(const vector<double>* ttags, const vector<double>* htags, const vector<double>* vtags, const vector<double>* offsets) { 
    const uint offset_len = offsets->size();
    vector<cc_point> result(offset_len);
    
    #pragma omp parallel for
    for(uint i=0; i<offset_len; ++i) {
        double os = offsets->at(i);
        
        cc_point ccdata;
        
        // H
        vector<double> tmptags = *htags;
        std::ranges::for_each(tmptags.begin(), tmptags.end(), [&os](double &ht) { ht= round((ht-os)*10)/10; });
        vector<double> cc_tags = {};
        cc_tags.reserve(4096);
        set_intersection(ttags->begin(), ttags->end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));

        ccdata.cc_h = cc_tags.size();
        ccdata.cc_h_tags = cc_tags;
        
        // V
        tmptags = *vtags;
        std::ranges::for_each(tmptags.begin(), tmptags.end(), [&os](double &vt) { vt = round((vt-os)*10)/10; });
        
        cc_tags = {};
        set_intersection(ttags->begin(), ttags->end(),
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
vector<cc_point> find_coincidences_fpga(const vector<double>* ttags, const vector<double>* htags, const vector<double>* vtags, const vector<double>* offsets, const vector<double>* ftags) {
    //round fpga tags to 10 ns
    vector<double> fpga_tags_rounded = *ftags;
    std::ranges::for_each(fpga_tags_rounded.begin(), fpga_tags_rounded.end(), [](double &ft) { ft= round(ft); });

    const uint offset_len = offsets->size();
    vector<cc_point> result(offset_len);

    #pragma omp parallel for
    for(uint i=0; i<offset_len; ++i) {
        double os = offsets->at(i);

        cc_point ccdata;

        // H
        vector<double> tmptags = *htags;
        std::ranges::for_each(tmptags.begin(), tmptags.end(), [&os](double &ht) { ht= round((ht-os)*10)/10; });
        vector<double> cc_tags = {};
        cc_tags.reserve(4096);
        set_intersection(ttags->begin(), ttags->end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));

        //check coincidences between fpga and cc_tags
        vector<double> ccc_tags = {};
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
        std::ranges::for_each(tmptags.begin(), tmptags.end(), [&os](double &vt) { vt = round((vt-os)*10)/10; });

        cc_tags = {};
        set_intersection(ttags->begin(), ttags->end(),
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
void cc_point_to_ccstruct(const vector<cc_point> *pts, ccstruct *ccs) {
    uint numpts = pts->size();
    vector<long> cc_h = {};
    cc_h.reserve(numpts);
    vector<long> cc_v = {};
    cc_v.reserve(numpts);
    vector<vector<double>> cc_h_tags = {{}};
    cc_h_tags.reserve(numpts);
    vector<vector<double>> cc_v_tags = {{}};
    cc_v_tags.reserve(numpts);
    vector<double> offsets = {};
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
*** create vector holding range of values
*/
template<typename T>
vector<T> arange(const T start, const T stop, const T step) {
    vector<T> values;
    for (T value = start; value < stop; value += step)
        values.emplace_back(value);
    return values;
}




/********************************************************************************
*** find and replace in a string (find 1, NOT N)
*/
bool stringreplace(string& str, const string& from, const string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

/********************************************************************************
*** print contents of a ccstruct
*/
void print_ccstruct(const ccstruct* dat) {
    cout << "offsets" << endl;
    for (uint i=0; i<dat->offsets.size(); ++i) {
        cout << dat->offsets[i] << ", ";
    }
    cout << endl;
    
    cout << "cc_h" << endl;
    for (uint i=0; i<dat->cc_h.size(); ++i) {
        cout << dat->cc_h[i] << ", ";
    }
    cout << endl;
    
    cout << "cc_v" << endl;
    for (uint i=0; i<dat->cc_v.size(); ++i) {
        cout << dat->cc_v[i] << ", ";
    }
    cout << endl;
    
    cout << "cc_h_tags" << endl;
    for (uint i=0; i<dat->cc_h_tags.size(); ++i) {
        for (uint j=0; j<dat->cc_h_tags[i].size(); ++j) {
            cout << dat->cc_h_tags[i][j] << ", ";
        }
        cout << ", ";
    }
    cout << endl;
    
    cout << "cc_v_tags" << endl;
    for (uint i=0; i<dat->cc_v_tags.size(); ++i) {
        for (uint j=0; j<dat->cc_v_tags[i].size(); ++j) {
            cout << dat->cc_v_tags[i][j] << ", ";
        }
         cout << ", ";
    }
    
    cout << "meastime: " << dat->meastime << endl;
    cout << "number of trigger tags: " << dat->num_trigger_tags << endl;
    cout << "number of fpga tags: " << dat->num_fpga_tags << endl;
    cout << endl;
}

/********************************************************************************
*** save a ccstruct to disk using protobuf
*/
int ccstruct_protobuf_todisk(const ccstruct* data, const string fname) {
    cout << "writing to file " << fname << endl;
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
        fstream ofs(fname, ios::out | ios::trunc | ios::binary);
        if (!pbdat.SerializeToOstream(&ofs)) {
            cerr << "Failed to write." << endl;
            return -1;
        }
        ofs.close();
    }
    
    return 0;
}


/********************************************************************************
*** load a ccstruct from disk using protobuf
*/
int ccstruct_protobuf_fromdisk(ccstruct* data, const string fname)
{
    //read from disk
    ccset::ccset_data read_pbdat;
    {
        fstream ifs(fname, ios::in | ios::binary);
        if (!read_pbdat.ParseFromIstream(&ifs)) {
            cerr << "Failed to parse." << endl;
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
        vector<double> tmparr = { };
        for (int j=0; j<read_pbdat.cc_h_tags(i).arr_size(); ++j) {
            tmparr.emplace_back(read_pbdat.cc_h_tags(i).arr(j));
        }
        data->cc_h_tags.emplace_back(tmparr);
    }
    for (int i=0; i<read_pbdat.cc_v_tags_size(); ++i) {
        vector<double> tmparr = { };
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
vector<string> get_new_tagfiles() {
    string path = std::filesystem::current_path();
    vector<string> alltagfiles = {};
    vector<string> analyzedfiles = {};
    vector<string> newtagfiles = {};
    
    string subtruncstr;
    if (cfg.TRUNCATE_S > 0) {
        string truncstr = string(3 - to_string(cfg.TRUNCATE_S).length(), '0') + to_string(cfg.TRUNCATE_S);
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
            if (p.path().filename().string().find(subtruncstr) != string::npos) {
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
            string h5string = *it;
            stringreplace(h5string, ".h5","");
            for ( auto p : analyzedfiles ) {
                string tmpstring = p;
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
*** check if a file exists
*/
bool fileExists(const std::string& fn)
{
    struct stat buf;
    if (stat(fn.c_str(), &buf) != -1) {
        return true;
    }
    return false;
}

/********************************************************************************
*** read config file
*/
void read_config() {
    ifstream cfgf("cc.cfg");
    if (cfgf.is_open()) {
        string line;
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
                cout << "unknown config name: " << name << endl;
            }
        }
        cfgf.close();
    } else {
        cerr << "Couldn't open config file for reading.\n";
    }

    cout << "---config---" << endl;
    cout << "hist_start    : " << cfg.HIST_START << "ns" << endl;
    cout << "hist_stop     : " << cfg.HIST_STOP  << "ns" << endl;
    cout << "hist_step     : " << cfg.HIST_STEP  << "ns" << endl;
    cout << "thigger chn   : " << cfg.CHN_TR     << endl;
    cout << "fpga chn      : " << cfg.CHN_FPGA   << endl;
    cout << "H chn         : " << cfg.CHN_H      << endl;
    cout << "V chn         : " << cfg.CHN_V      << endl;
    cout << "truncate at   : " << cfg.TRUNCATE_S << "s" << endl;
    cout << "use fpga chan : " << cfg.FPGA_USE   << endl;
    cout << "fpga delay min: " << cfg.FPGA_DELAY_MIN << "ns" << endl;
    cout << "fpga delay max: " << cfg.FPGA_DELAY_MAX << "ns" << endl;
    cout << "patterns:     : " << endl;
    for (auto e: cfg.patterns) {
        print_vector(e);
    }

}

/********************************************************************************
*** parse config string containing patternspatterns
*/
std::vector<std::vector<uint16_t>> parse_patterns(const std::string instring) {
    std::string s = instring;
    if (!(s.starts_with("{{"))) {
        cout << "pattern string malformatted: " << s << endl;
    }
    if (!(s.ends_with("}}"))) {
        cout << "pattern string malformatted: " << s << endl;
    }
    
    s.erase(0,1);
    s.erase(s.length()-1,1);
    
    size_t vec_start_idx = 0;
    size_t vec_stop_idx = 0;
    std::vector<size_t> sep_idxs = {};
    size_t curr_sep_idx = 0;
    
    std::string vector_string = "";
    std::vector<uint16_t> pattern = {};
    std::vector<std::vector<uint16_t>> patterns = {};
    std::string channelstr = "";
    // separate patterns
    while (vec_stop_idx < s.length()-1) {
        vec_start_idx =  s.find("{", vec_start_idx);
        vec_stop_idx = s.find("}", vec_start_idx);
        vector_string = s.substr(vec_start_idx, vec_stop_idx-vec_start_idx+1);
        
        sep_idxs = {};
        sep_idxs.push_back(vector_string.find('{'));
        curr_sep_idx = vector_string.find(',', 0);
        while(curr_sep_idx != string::npos) {
            sep_idxs.push_back(curr_sep_idx);
            curr_sep_idx = vector_string.find(',',curr_sep_idx+1);
        }
        sep_idxs.push_back(vector_string.find('}'));
        
        //get channels
        pattern = {};
        for (size_t i=0; i< sep_idxs.size()-1; ++i) {
            channelstr = vector_string.substr(sep_idxs[i]+1, sep_idxs[i+1]-sep_idxs[i]-1);
            channelstr = std::regex_replace(channelstr, std::regex(R"([\D])"), "");
            pattern.push_back(stoi(channelstr));
        }
        patterns.push_back(pattern);
        
        vec_start_idx = vec_stop_idx;
    }
    
    return patterns;
}

/********************************************************************************
*** create vector holding range of values
*/
template<typename T>
void print_vector(const std::vector<T> v) {
    for (T e : v) {
        if (e != v.back()) {
            std::cout << e << ", ";
        } else {
            std::cout << e << endl;
        }
    }
}

/********************************************************************************
*** turn vector of vectors into vector
*/
template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> vv) {
    std::vector<T> result;
    for (auto v: vv) {
        result.insert(result.end(), v.begin(), v.end());
    }
    return result;
}

/********************************************************************************
*** find unique entries in vector
*/
template<typename T>
std::vector<T> unique(const std::vector<T> v) {
    std::vector<T> result;
    for (auto e: v) {
        if (std::find(result.begin(), result.end(), e) != result.end()) {
            result.push_back(e);
        }
    }
    return result;
}
