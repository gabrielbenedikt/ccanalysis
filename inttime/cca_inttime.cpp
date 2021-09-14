#include "cca_inttime.h"

int main(void)
{   
    read_config();
    
    vector<string> newtagfiles = get_new_tagfiles();
    
    cout << "new files: "<< endl;
    for (auto e: newtagfiles) {
        cout << e << endl;
    }
    
//     #pragma omp parallel for
    for ( uint i = 0; i<newtagfiles.size(); ++i ) {
        auto starttime = high_resolution_clock::now();
        string fn = newtagfiles[i];
        string datasetPath = "/tags/block0_values";
        
        /*
        * read hdf file
        */
        long data_len = 0;
        long *data = readHDFfile(fn, datasetPath, data_len);
        
        /*
        * separate into channels
        * in internal clock time
        */
        vector<long> cnt_t = {};
        vector<long> cnt_h = {};
        vector<long> cnt_v = {};
        separate_tags_per_channels(data, data_len, cnt_t, cnt_h, cnt_v);
        
        
        //round trigger tags
        for(auto& tt : cnt_t) {
            tt=roundUp(tt, BIN_WIDTH);
        }
        
        cout << "file " << i+1 << " of " << newtagfiles.size() << ": filename: " << fn << endl;
        cout << "number of trigger / H / V tags: " << cnt_t.size() << " / " << cnt_h.size() << " / " << cnt_v.size() << endl;
        
        /*
        * find coincidences
        */
        double meastime = (cnt_t.back()-cnt_t[0])*pow(10,-9);
        cout << "measurement time: " << meastime << endl;
        vector<long> offsets = arange<long>(BIN_START,BIN_STOP,BIN_WIDTH);
        
        vector<cc_point> hist_points = find_coincidences(&cnt_t, &cnt_h, &cnt_v, &offsets);
        
        //convert point vector to single struct
        ccstruct hist_combined;
        cc_point_to_ccstruct(&hist_points, &hist_combined);
        hist_combined.meastime = meastime;
        
        // save analysis to disk
        string savefname = fn;
        stringreplace(savefname,".h5","_ccs.pbdat");
        ccstruct_protobuf_todisk(&hist_combined, savefname);
        
        auto endtime = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(endtime - starttime).count();
        
        SPENT_TIME += duration; 
        ETA = (SPENT_TIME/i)*(newtagfiles.size()-i);
        cout << "Analysis time / average analysis time / ETA: " << duration << "s / " << SPENT_TIME/i <<"s / " << ETA/60 <<"m" << endl;
        
        delete[] data;
    }
    
    return 0;
}

/********************************************************************************
*** read hdf file
*/
long* readHDFfile(string fn, string datasetPath, long& out_data_len)
{
    
    // Open HDF5 file handle, read only
    H5File file(fn.c_str(), H5F_ACC_RDONLY);

    // access the required dataset by path name
    DataSet dataset = file.openDataSet(datasetPath.c_str());

    // get the dataspace
    DataSpace dataspace = dataset.getSpace();
  
    
//    H5T_class_t type_class = dataset.getTypeClass();
//   /*
//    * Get class of datatype and print message if it's an integer.
//    */
//    if( type_class == H5T_INTEGER ) {
//        cout << "Data set has INTEGER type" << endl;
//        /*
//        * Get the integer datatype
//        */
//        IntType intype = dataset.getIntType();
//        /*
//        * Get order of datatype and print message if it's a little endian.
//        */
//        H5std_string order_string;
//        H5T_order_t order = intype.getOrder( order_string );
//        cout << order_string << endl;
//        /*
//        * Get size of the data element stored in file and print it.
//        */
//       size_t size = intype.getSize();
//       cout << "Data size is " << size << endl;
//   }
    
    /*
    * Get the number of dimensions in the dataspace.
    */
    //int rank = dataspace.getSimpleExtentNdims();
    /*
    * Get the dimension size of each dimension in the dataspace and
    * display them.
    */
    hsize_t dims_out[2];
    (void)(dataspace.getSimpleExtentDims( dims_out, NULL));
    //cout << "rank " << rank << ", dimensions " <<
    //(unsigned long)(dims_out[0]) << " x " <<
    //(unsigned long)(dims_out[1]) << endl;

    int NX=dims_out[0];
    int NY=dims_out[1];
    /*
    * Output buffer initialization.
    */
    long* data_out = new long[NX*NY];
    
    out_data_len = NX*NY;
    /*
    * Read data from hyperslab in the file into the hyperslab in
    * memory and display the data.
    */
    dataset.read( data_out, PredType::NATIVE_LONG);
    
    // close the HDF5 file
    file.close();
    
    return data_out;
}


/********************************************************************************
*** separate tag vectors
*/
void separate_tags_per_channels(long* tags, long numtags, vector<long>& out_countst, vector<long>& out_countsh, vector<long>& out_countsv)
{   
    vector<long> cnt_t = {  };
    vector<long> cnt_h = {  };
    vector<long> cnt_v = {  };
    
    for (int i=0; i<numtags; i+=2) {
        if (tags[i]==CHN_T) {
            cnt_t.push_back(tags[i+1]);
        } else if (tags[i]==CHN_H) {
            cnt_h.push_back(tags[i+1]);
        } else if (tags[i]==CHN_V) {
            cnt_v.push_back(tags[i+1]);
        }
    }
    
    out_countst = cnt_t;
    out_countsh = cnt_h;
    out_countsv = cnt_v;
}


/********************************************************************************
*** find coincidences - STL BST
*/
vector<cc_point> find_coincidences(vector<long>* ttags, vector<long>* htags, vector<long>* vtags, vector<long>* offsets) {
    vector<cc_point> result(offsets->size());
    
    #pragma omp parallel for
    for(uint i=0; i<offsets->size(); ++i) {
        long os = offsets->at(i);
        
        vector<long> tmptagsh = *htags;
        vector<long> tmptagsv = *vtags;
        int cc_h;
        int cc_v;
        
        for(auto& ht : tmptagsh) {
            roundUp(ht-os, BIN_WIDTH);
        }
        for(auto& vt : tmptagsv) {
            roundUp(vt-os, BIN_WIDTH);
        }
        
        vector<long> cc_h_tags = {};
        vector<long> cc_v_tags = {};
        
        for (long htag : tmptagsh) {
            if (binary_search(ttags->begin(), ttags->end(), htag)) {
                cc_h_tags.push_back(htag);
            }
        }
        for (long vtag : tmptagsv) {
            if (binary_search(ttags->begin(), ttags->end(), vtag)) {
                cc_v_tags.push_back(vtag);
            }
        }
                    
        cc_h = cc_h_tags.size();
        cc_v = cc_v_tags.size();
        
        cc_point ccdata;
        ccdata.cc_h = cc_h;
        ccdata.cc_v = cc_v;
        ccdata.cc_h_tags = cc_h_tags;
        ccdata.cc_v_tags = cc_v_tags;
        ccdata.offset = os;
        
        result[i] = ccdata;
    }
    
    return result;
}

//********************************************************************************
//*** find coincidences - custom binary search
//*/
vector<cc_point> find_coincidences2(vector<long>* ttags, vector<long>* htags, vector<long>* vtags, vector<long>* offsets) {
    vector<cc_point> result(offsets->size());
    
    #pragma omp parallel for
    for(uint i=0; i<offsets->size(); ++i) {
        long os = offsets->at(i);
        
        vector<long> tmptagsh = *htags;
        vector<long> tmptagsv = *vtags;
        int cc_h;
        int cc_v;
        
        for(auto& ht : tmptagsh) {
            roundUp(ht-os, BIN_WIDTH);
        }
        for(auto& vt : tmptagsv) {
            roundUp(vt-os, BIN_WIDTH);
        }
        
        vector<long> cc_h_tags = {};
        vector<long> cc_v_tags = {};
        
        int r = 0;
        int l = 0;
        int mid = 0;
        int minidx=0;
        for (long htag: tmptagsh) {
            l = minidx;
            r = ttags->size()-1;
            
            while (l<=r) {
                mid = l + floor((r - l)/2);
                if (ttags->at(mid) == htag) {
                    cc_h_tags.push_back(htag);
                    minidx=mid;
                    break;
                } else if (ttags->at(mid) < htag) {
                    l = mid + 1;
                } else {
                    r = mid - 1;
                }
            }
        }
                
        r = 0;
        l = 0;
        mid = 0;
        minidx=0;
        for (long vtag: tmptagsv) {
            l = minidx;
            r = ttags->size()-1;
            
            while (l<=r) {
                mid = l + floor((r - l)/2);
                if (ttags->at(mid) == vtag) {
                    cc_v_tags.push_back(vtag);
                    minidx=mid;
                    break;
                } else if (ttags->at(mid) < vtag) {
                    l = mid + 1;
                } else {
                    r = mid - 1;
                }
            }
        }
                    
        cc_h = cc_h_tags.size();
        cc_v = cc_v_tags.size();
        
        cc_point ccdata;
        ccdata.cc_h = cc_h;
        ccdata.cc_v = cc_v;
        ccdata.cc_h_tags = cc_h_tags;
        ccdata.cc_v_tags = cc_v_tags;
        ccdata.offset = os;
        
        result[i] = ccdata;
       
    }
    
    return result;
}


//********************************************************************************
//*** find coincidences - unordered set
//*/
vector<cc_point> find_coincidences3(vector<long>* ttags, vector<long>* htags, vector<long>* vtags, vector<long>* offsets) {
    vector<cc_point> result(offsets->size());
    
    std::unordered_set<long> tmpttags(ttags->begin(), ttags->end());
    #pragma omp parallel for
    for(uint i=0; i<offsets->size(); ++i) {
        long os = offsets->at(i);
        
        vector<long> tmptagsh = *htags;
        vector<long> tmptagsv = *vtags;
        int cc_h;
        int cc_v;
        
        for(auto& ht : tmptagsh) {
            roundUp(ht-os, BIN_WIDTH);
        }
        for(auto& vt : tmptagsv) {
            roundUp(vt-os, BIN_WIDTH);
        }
        
        vector<long> cc_h_tags = {};
        vector<long> cc_v_tags = {};
        
        for (long htag : tmptagsh) {
            if (tmpttags.contains(htag)) {
                cc_h_tags.push_back(htag);
            }
        }
        for (long vtag : tmptagsv) {
            if (tmpttags.contains(vtag)) {
                cc_v_tags.push_back(vtag);
            }
        }
                    
        cc_h = cc_h_tags.size();
        cc_v = cc_v_tags.size();
        
        cc_point ccdata;
        ccdata.cc_h = cc_h;
        ccdata.cc_v = cc_v;
        ccdata.cc_h_tags = cc_h_tags;
        ccdata.cc_v_tags = cc_v_tags;
        ccdata.offset = os;
        
        result[i] = ccdata;
    }
    
    return result;
}

//********************************************************************************
//*** find coincidences - linear search
//*/
vector<cc_point> find_coincidences4(vector<long>* ttags, vector<long>* htags, vector<long>* vtags, vector<long>* offsets) {
    vector<cc_point> result(offsets->size());
    
    #pragma omp parallel for
    for(uint i=0; i<offsets->size(); ++i) {
        long os = offsets->at(i);
        
        vector<long> tmptagsh = *htags;
        vector<long> tmptagsv = *vtags;
        int cc_h;
        int cc_v;
        
        for(auto& ht : tmptagsh) {
            roundUp(ht-os, BIN_WIDTH);
        }
        for(auto& vt : tmptagsv) {
            roundUp(vt-os, BIN_WIDTH);
        }
        
        vector<long> cc_h_tags = {};
        vector<long> cc_v_tags = {};
        
        uint startidx=0;
        uint endidx=ttags->size();
        for (long htag : tmptagsh) {
            for (uint idx = startidx; i<endidx; ++i ) {
                if (htag == ttags->at(i)) {
                    cc_h_tags.push_back(htag);
                    startidx = idx;
                    break;
                } else if (htag < ttags->at(i)) {
                    startidx = idx;
                    break;
                }
            }
        }
        
        startidx=0;
        for (long vtag : tmptagsv) {
            for (uint idx = startidx; i<endidx; ++i ) {
                if (vtag == ttags->at(i)) {
                    cc_v_tags.push_back(vtag);
                    startidx = idx;
                    break;
                } else if (vtag < ttags->at(i)) {
                    startidx = idx;
                    break;
                }
            }
        }
                    
        cc_h = cc_h_tags.size();
        cc_v = cc_v_tags.size();
        
        cc_point ccdata;
        ccdata.cc_h = cc_h;
        ccdata.cc_v = cc_v;
        ccdata.cc_h_tags = cc_h_tags;
        ccdata.cc_v_tags = cc_v_tags;
        ccdata.offset = os;
        
        result[i] = ccdata;
    }
    
    return result;
}

/********************************************************************************
*** separate tag vectors
*/
void cc_point_to_ccstruct(vector<cc_point> *pts, ccstruct *ccs) {
    vector<long> cc_h = {};
    vector<long> cc_v = {};
    vector<vector<long>> cc_h_tags = {{}};
    vector<vector<long>> cc_v_tags = {{}};
    vector<long> offsets = {};
    for (auto pnt: *pts) {
        cc_h.push_back(pnt.cc_h);
        cc_v.push_back(pnt.cc_v);
        cc_h_tags.push_back(pnt.cc_h_tags);
        cc_v_tags.push_back(pnt.cc_v_tags);
        offsets.push_back(pnt.offset);
    }
    
    ccs->offsets = offsets;
    ccs->cc_h = cc_h;
    ccs->cc_v = cc_v;
    ccs->cc_h_tags = cc_h_tags;
    ccs->cc_v_tags = cc_v_tags;
}


/********************************************************************************
*** separate tag vectors
*/
template<typename T>
vector<T> arange(T start, T stop, T step) {
    vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
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
void print_ccstruct(ccstruct* dat) {
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
    cout << endl;
}

/********************************************************************************
*** save a ccstruct to disk using protobuf
*/
int ccstruct_protobuf_todisk(ccstruct* data, string fname) {
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
int ccstruct_protobuf_fromdisk(ccstruct* data, string fname) 
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
        data->offsets.push_back(read_pbdat.offsets(i));
    }
    for (int i=0; i<read_pbdat.cc_h_size(); ++i) {
        data->cc_h.push_back(read_pbdat.cc_h(i));
    }
    for (int i=0; i<read_pbdat.cc_v_size(); ++i) {
        data->cc_v.push_back(read_pbdat.cc_v(i));
    }
    for (int i=0; i<read_pbdat.cc_h_tags_size(); ++i) {
        vector<long> tmparr = { };
        for (int j=0; j<read_pbdat.cc_h_tags(i).arr_size(); ++j) {
            tmparr.push_back(read_pbdat.cc_h_tags(i).arr(j));
        }
        data->cc_h_tags.push_back(tmparr);
    }
    for (int i=0; i<read_pbdat.cc_v_tags_size(); ++i) {
        vector<long> tmparr = { };
        for (int j=0; j<read_pbdat.cc_v_tags(i).arr_size(); ++j) {
            tmparr.push_back(read_pbdat.cc_v_tags(i).arr(j));
        }
        data->cc_v_tags.push_back(tmparr);
    }
    data->meastime=read_pbdat.meastime();
    
    return 0;
}

/********************************************************************************
*** find new tagfiles
*/
vector<string> get_new_tagfiles() {
    string path = fs::current_path();
    vector<string> alltagfiles = {};
    vector<string> analyzedfiles = {};
    vector<string> newtagfiles = {};
    for ( auto & p : fs::directory_iterator(path) ) {
        if ( p.path().extension() == ".h5" ) {
            if (p.path().filename() != "tomos.h5") {
                alltagfiles.push_back(p.path());
            }
        } else if ( p.path().extension() == ".pbdat" ) {
            analyzedfiles.push_back(p.path());
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
                stringreplace(tmpstring,"_ccs.pbdat","");
                
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
    ifstream cfgf("ccint.cfg");
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
            
            if (name == "bin_start") {
                BIN_START = stoi(value);
            } else if (name == "bin_stop") {
                BIN_STOP = stoi(value);
            } else if (name == "bin_width") {
                BIN_WIDTH = stoi(value);
            } else if (name == "chn_h") {
                CHN_H = stoi(value);
            } else if (name == "chn_v") {
                CHN_V = stoi(value);
            } else if (name == "chn_t") {
                CHN_T = stoi(value);
            } else {
                cout << "unknown config name: " << name << endl;
            }
        }
        cfgf.close();
    } else {
        cerr << "Couldn't open config file for reading.\n";
    }
    
    cout << "---config---" << endl;
    cout << "bin_start  : " << BIN_START  << endl;
    cout << "bin_stop   : " << BIN_STOP   << endl;
    cout << "bin_width  : " << BIN_WIDTH  << endl;
    cout << "thigger chn: " << CHN_T      << endl;
    cout << "H chn      : " << CHN_H      << endl;
    cout << "V chn      : " << CHN_V      << endl;
}


/********************************************************************************
*** rounds number to next hight multiple of arbitrary number
*/
long roundUp(long numToRound, long multiple)
{
    if (multiple == 0)
        return numToRound;

    long remainder = abs(numToRound) % multiple;
    if (remainder == 0)
        return numToRound;

    if (numToRound < 0)
        return -(abs(numToRound) - remainder);
    else
        return numToRound + multiple - remainder;
}
