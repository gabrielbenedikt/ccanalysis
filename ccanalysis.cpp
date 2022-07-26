#include "ccanalysis.h"

using namespace std;
//using namespace H5;

int main(void)
{   
    omp_set_num_threads(8);

    read_config();
    
    vector<string> newtagfiles = get_new_tagfiles();
    
    cout << "new files: "<< endl;
    for (const auto &e: newtagfiles) {
        cout << e << endl;
    }
    
    uint nanalyzed = 0;
    //#pragma omp parallel for num_threads(2)
    for ( uint i = 0; i<newtagfiles.size(); ++i ) {
        auto starttime = std::chrono::high_resolution_clock::now();
        string fn = newtagfiles[i];
        string datasetPath = "/tags/block0_values";
        cout << "file " << nanalyzed << " of " << newtagfiles.size() << ": filename: " << fn << endl;
        
        /*
        * read hdf file
        */
        long long data_len = 0;
        long long *data = readHDFfile(fn, datasetPath, data_len);
        
        /*
        * separate into channels
        * already converted to ns
        */
        vector<double> cnt_tr = {};
        vector<double> cnt_h = {};
        vector<double> cnt_v = {};
        vector<double> cnt_fpga = {};
        separate_tags_per_channels(data, data_len, cnt_tr, cnt_h, cnt_v, cnt_fpga);
        
        std::for_each(std::execution::par_unseq, cnt_tr.begin(), cnt_tr.end(), [](double &tt) { tt=std::round(tt*10)/10; });
#if defined(PREPROCESS_TRIGGER_WITH_FPGA_TAGS) || defined(CC_ALGO_STL_SETINTERSECT_TRIGGER)
        std::for_each(std::execution::par_unseq, cnt_fpga.begin(), cnt_fpga.end(), [](double &ft) { ft=std::round(ft*10)/10; });
#endif

        
#ifdef PREPROCESS_TRIGGER_WITH_FPGA_TAGS
        cnt_tr = find_coincidences_trigger_fpga(&cnt_fpga, &cnt_tr);
#endif
        /*
        * find coincidences
        */
        double meastime = (cnt_tr.back()-cnt_tr.front())*pow(10,-9);
        vector<double> offsets = arange<double>(HIST_START,HIST_STOP,HIST_STEP);

#ifdef CC_ALGO_STL_SETINTERSECT_TRIGGER
            vector<cc_point> hist_points = find_coincidences(&cnt_tr, &cnt_h, &cnt_v, &offsets, &cnt_fpga);
#else
            vector<cc_point> hist_points = find_coincidences(&cnt_tr, &cnt_h, &cnt_v, &offsets);
#endif
            
        //convert point vector to single struct
        ccstruct hist_combined;
        cc_point_to_ccstruct(&hist_points, &hist_combined);
        hist_combined.meastime = meastime;
        hist_combined.num_trigger_tags = cnt_tr.size();
        hist_combined.num_fpga_tags = cnt_fpga.size();
        
        // save analysis to disk
        string savefname = fn;
        string newext;
        if (TRUNCATE_S > 0) {
            string truncstr = string(3 - to_string(TRUNCATE_S).length(), '0') + to_string(TRUNCATE_S);
            newext = "_trunc"+ truncstr +"_ccs.pbdat";
        } else {
            newext = "_trunc000_ccs.pbdat";
        }
        stringreplace(savefname,".h5", newext);
        ccstruct_protobuf_todisk(&hist_combined, savefname);
        
        auto endtime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endtime - starttime).count();
        ++nanalyzed;
        SPENT_TIME += duration; 
        ETA = (SPENT_TIME/nanalyzed)*(newtagfiles.size()-nanalyzed);
        
        cout << "number of trigger / fpga / H / V tags: " << cnt_fpga.size() << " / " << cnt_tr.size() << " / " << cnt_h.size() << " / " << cnt_v.size() << endl;
        cout << "measurement time: " << meastime << endl;
        cout << "analysis T / avg T / ETA: " << duration << "s / " << SPENT_TIME/i <<"s / " << ETA/(60*omp_get_num_threads()) <<"m" << endl;
        
        delete[] data;
    }
    
    cout << "done" << endl;
    
    return 0;
}

/********************************************************************************
*** read hdf file
*/
long long* readHDFfile(const string fn, const string datasetPath, long long& out_data_len)
{
    
    // Open HDF5 file handle, read only
    H5::H5File file(fn.c_str(), H5F_ACC_RDONLY);

    // access the required dataset by path name
    H5::DataSet dataset = file.openDataSet(datasetPath.c_str());

    // get the dataspace
    H5::DataSpace dataspace = dataset.getSpace();
  
    
//    H5T_class_t type_class = dataset.getTypeClass();
   /*
    * Get class of datatype and print message if it's an integer.
    */
//    if( type_class == H5T_INTEGER ) {
//        cout << "Data set has INTEGER type" << endl;
//        H5::IntType intype = dataset.getIntType(); // Get the integer datatype
//
//        H5std_string order_string; // Get order of datatype and print endianness.
//        /*H5T_order_t order = */ intype.getOrder( order_string );
//        cout << order_string << endl;
//
//       size_t size = intype.getSize(); //Get size of the data element stored in file and print it.
//       cout << "Data size is " << size << endl;
//
//    } else if ( type_class == H5T_FLOAT ) {
//        cout << "Data set has FLOAT type" << endl;
//    } else {
//        cout << "Data set type: " << type_class << endl;
//    }
    
    /*
    * Get the number of dimensions in the dataspace.
    */
//    int rank = dataspace.getSimpleExtentNdims();
    /*
    * Get the dimension size of each dimension in the dataspace and
    * display them.
    */
    hsize_t dims_out[2];
    /*int dimensions = */ dataspace.getSimpleExtentDims( dims_out, NULL);
//    cout << "rank " << rank << ", dimensions " <<
//    (unsigned long)(dims_out[0]) << " x " <<
//    (unsigned long)(dims_out[1]) << endl;

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
    //cout << "separate" << endl;
    if (TRUNCATE_S > 0) {
        //cout << "trunc" << endl;
        double firsttag = tags[1]*CS;
        double lasttag = 0;
        //cout << "first tag  : " << firsttag*pow(10,-9) << endl;
        //cout << "end tag    : " << tags[numtags-1]*CS*pow(10,-9) << endl;
        //cout << "trunkate at: " << TRUNCATE_S << endl;
        for (long long i=2; i<numtags; i+=2) {
            lasttag = tags[i+1]*CS;
            //some tags are zero. this breaks the comparison
            if (lasttag!=0) {
                if ((lasttag-firsttag)>TRUNCATE_S*pow(10,9)) {
                    //cout << "found tag to truncate at idx " << i << ": tag " << lasttag*pow(10,-9) << endl;
                    maxtag = i;
                    break;
                }
            }
        }
    }


    double currtag;
    for (long long i=0; i<maxtag; i+=2) {
        currtag = tags[i+1]*CS;
        //some tags are zero. this breaks stuff. exclude them. effing tagger.
        if (currtag!=0) {
            if (tags[i]==CHN_TR) {
                cnt_tr.emplace_back(tags[i+1]*CS);
            } else if (tags[i]==CHN_H) {
                cnt_h.emplace_back(tags[i+1]*CS);
            } else if (tags[i]==CHN_V) {
                cnt_v.emplace_back(tags[i+1]*CS);
            } else if (tags[i]==CHN_FPGA) {
                cnt_fpga.emplace_back(tags[i+1]*CS);
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
        //auto it = std::ranges::find_if(cnt_tr, is_neg);
        auto it = std::find_if(std::execution::par_unseq, cnt_tr.begin(), cnt_tr.end(), is_neg);
        
        // add the last positive element to all remaining elements in vector
        // by taking this difference, the first tag after the jump will be a 'double click'.
        // Better: find out range of tags and add this
        
        double diff = pow(2,42);                                    /////////////////////double diff = *prev(it) - *it;
        std::transform(std::execution::par_unseq, it, cnt_tr.end(), it, [&diff](double x) {return x+diff; } ); 
        
        //now, do the same for cnt_h and cnt_v
        if ((cnt_h.front() > 0) && (cnt_h.back() < 0)) {
            //auto it = std::ranges::find_if(cnt_h, is_neg);
            auto it = std::find_if(std::execution::par_unseq, cnt_h.begin(), cnt_h.end(), is_neg);
            //diff = *prev(it) - *it;
            std::transform(std::execution::par_unseq, it, cnt_h.end(), it, [&diff](double x) {return x+diff; } ); 
        }
        if ((cnt_v.front() > 0) && (cnt_v.back() < 0)) {
            //auto it = std::ranges::find_if(cnt_v, is_neg);
            auto it = std::find_if(std::execution::par_unseq, cnt_v.begin(), cnt_v.end(), is_neg);
            //diff = *prev(it) - *it;
            std::transform(std::execution::par_unseq, it, cnt_v.end(), it, [&diff](double x) {return x+diff; } ); 
        }
        if ((cnt_fpga.front() > 0) && (cnt_fpga.back() < 0)) {
            //auto it = std::ranges::find_if(cnt_fpga, is_neg);
            auto it = std::find_if(std::execution::par_unseq, cnt_fpga.begin(), cnt_fpga.end(), is_neg);
            //diff = *prev(it) - *it;
            std::transform(std::execution::par_unseq, it, cnt_fpga.end(), it, [&diff](double x) {return x+diff; } ); 
        }
    } // else same sign, thus no tag jump
}



/********************************************************************************
*** find coincidences - between trigger and fpga signal. return trimmed trigger tag vector
*/
std::vector<double> find_coincidences_trigger_fpga(const std::vector<double>* ftags, const std::vector<double>* ttags) {

    std::vector<double> trimmed_ttags;
    
    std::cout << "looking for trigger/fpga coincidences" << std::endl;
    
    double tt = 0;
    double ft=0;
    size_t startj = 0;
    
    for (size_t i = 0; i < ttags->size(); ++i) {
        tt = ttags->at(i);
        for (size_t j = startj; j < ftags->size(); ++j) {
            ft = ftags->at(j);
            
            if (ft > (tt+FPGA_T_MIN)) {
                if (ft < (tt+FPGA_T_MAX)) {
                    trimmed_ttags.emplace_back(tt);
                    startj = j;
                }
                break;
            } else {
                startj = j;
            }
        }
    }
    
    std::cout << "done" << std::endl;
    std::cout << "#all  triggers: " << ttags->size() << std::endl;
    std::cout << "#fpga+triggers: " << trimmed_ttags.size() << std::endl;
    
    return trimmed_ttags;
}

#if CCALGO==STL_SETINTERSECT
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
        std::for_each(std::execution::unseq, tmptags.begin(), tmptags.end(), [&os](double &ht) { ht = std::round((ht-os)*10)/10; });
        vector<double> cc_tags = {};
        cc_tags.reserve(4096);
        std::set_intersection(std::execution::unseq, ttags->begin(), ttags->end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));

        ccdata.cc_h = cc_tags.size();
        ccdata.cc_h_tags = cc_tags;
        
        // V
        tmptags = *vtags;
        std::for_each(std::execution::unseq, tmptags.begin(), tmptags.end(), [&os](double &vt) { vt = std::round((vt-os)*10)/10; });
        
        cc_tags = {};
        std::set_intersection(std::execution::unseq, ttags->begin(), ttags->end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));

        ccdata.cc_v = cc_tags.size();
        ccdata.cc_v_tags = cc_tags;
        ccdata.offset = os;
        result[i] = ccdata;
    }
    return result;         
}

#elif CCALGO==STL_SETINTERSECT_TRIGGER
/********************************************************************************
*** find coincidences with trigger - set inersect - keep only coincidences where fpga signal is present too
*/
vector<cc_point> find_coincidences(const vector<double>* ttags, const vector<double>* htags, const vector<double>* vtags, const vector<double>* offsets, const vector<double>* ftags) {
    const uint offset_len = offsets->size();
    vector<cc_point> result(offset_len);
    
    #pragma omp parallel for
    for(uint i=0; i<offset_len; ++i) {
        double os = offsets->at(i);
        
        cc_point ccdata;
        
        // H
        vector<double> tmptags = *htags;
        
        std::for_each(std::execution::unseq, tmptags.begin(), tmptags.end(), [&os](double &ht) { ht= std::round((ht-os)*10)/10; });
        vector<double> cc_tags = {};
        cc_tags.reserve(4096);
        std::set_intersection(std::execution::unseq, ttags->begin(), ttags->end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));

        //check if trigger is present
        //TODO
        
        ccdata.cc_h = cc_tags.size();
        ccdata.cc_h_tags = cc_tags;
        
        // V
        tmptags = *vtags;
        std::for_each(std::execution::unseq, tmptags.begin(), tmptags.end(), [&os](double &vt) { vt = std::round((vt-os)*10)/10; });
        
        cc_tags = {};
        std::set_intersection(std::execution::unseq, ttags->begin(), ttags->end(),
                         tmptags.begin(), tmptags.end(),
                         std::back_inserter(cc_tags));

        ccdata.cc_v = cc_tags.size();
        ccdata.cc_v_tags = cc_tags;
        ccdata.offset = os;
        result[i] = ccdata;
    }
    return result;         
}

#elif CCALGO==BST
//********************************************************************************
//*** find coincidences - custom binary search
//*/
vector<cc_point> find_coincidences(const vector<double>* ttags, const vector<double>* htags, const vector<double>* vtags, const vector<double>* offsets) {
    vector<cc_point> result(offsets->size());
    
    #pragma omp parallel for
    for(uint i=0; i<offsets->size(); ++i) {
        double os = offsets->at(i);
        
        vector<double> tmptagsh = *htags;
        vector<double> tmptagsv = *vtags;
        int cc_h;
        int cc_v;
        std::ranges::for_each(tmptagsh.begin(), tmptagsh.end(), [&os](double &ht) { ht= round((ht-os)*10)/10; });
        std::ranges::for_each(tmptagsv.begin(), tmptagsv.end(), [&os](double &vt) { vt= round((ht-os)*10)/10; });
        
        vector<double> cc_h_tags = {};
        vector<double> cc_v_tags = {};
        
        int r = 0;
        int l = 0;
        int mid = 0;
        int minidx=0;
        for (double htag: tmptagsh) {
            l = minidx;
            r = ttags->size()-1;
            
            while (l<=r) {
                mid = l + floor((r - l)/2);
                if (ttags->at(mid) == htag) {
                    cc_h_tags.emplace_back(htag);
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
        for (double vtag: tmptagsv) {
            l = minidx;
            r = ttags->size()-1;
            
            while (l<=r) {
                mid = l + floor((r - l)/2);
                if (ttags->at(mid) == vtag) {
                    cc_v_tags.emplace_back(vtag);
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

#elif CCALGO==UNORDERED_SET
//********************************************************************************
//*** find coincidences - unordered set
//*/
vector<cc_point> find_coincidences(const vector<double>* ttags, const vector<double>* htags, const vector<double>* vtags, const vector<double>* offsets) {
    vector<cc_point> result(offsets->size());
    
    std::unordered_set<double> tmpttags(ttags->begin(), ttags->end());
    #pragma omp parallel for
    for(uint i=0; i<offsets->size(); ++i) {
        double os = offsets->at(i);
        
        vector<double> tmptagsh = *htags;
        vector<double> tmptagsv = *vtags;
        int cc_h;
        int cc_v;
        
        std::ranges::for_each(tmptagsh.begin(), tmptagsh.end(), [&os](double &ht) { ht= round((ht-os)*10)/10; });
        std::ranges::for_each(tmptagsv.begin(), tmptagsv.end(), [&os](double &vt) { vt= round((vt-os)*10)/10; });
        
        vector<double> cc_h_tags = {};
        vector<double> cc_v_tags = {};
        
        for (double htag : tmptagsh) {
            if (tmpttags.contains(htag)) {
                cc_h_tags.emplace_back(htag);
            }
        }
        for (double vtag : tmptagsv) {
            if (tmpttags.contains(vtag)) {
                cc_v_tags.emplace_back(vtag);
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

#elif CCALGO==LINEAR_SEARCH
//********************************************************************************
//*** find coincidences - linear search
//*/
vector<cc_point> find_coincidences(const vector<double>* ttags, const vector<double>* htags, const vector<double>* vtags, const vector<double>* offsets) {
    vector<cc_point> result(offsets->size());
    
    #pragma omp parallel for
    for(uint i=0; i<offsets->size(); ++i) {
        double os = offsets->at(i);
        
        vector<double> tmptagsh = *htags;
        vector<double> tmptagsv = *vtags;
        int cc_h;
        int cc_v;
        
        std::ranges::for_each(tmptagsh.begin(), tmptagsh.end(), [&os](double &ht) { ht= round((ht-os)*10)/10; });
        std::ranges::for_each(tmptagsv.begin(), tmptagsv.end(), [&os](double &vt) { vt= round((vt-os)*10)/10; });
        
        vector<double> cc_h_tags = {};
        vector<double> cc_v_tags = {};
        
        uint startidx=0;
        uint endidx=ttags->size();
        for (double htag : tmptagsh) {
            for (uint idx = startidx; i<endidx; ++i ) {
                if (htag == ttags->at(i)) {
                    cc_h_tags.emplace_back(htag);
                    startidx = idx;
                    break;
                } else if (htag < ttags->at(i)) {
                    startidx = idx;
                    break;
                }
            }
        }
        
        startidx=0;
        for (double vtag : tmptagsv) {
            for (uint idx = startidx; i<endidx; ++i ) {
                if (vtag == ttags->at(i)) {
                    cc_v_tags.emplace_back(vtag);
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

#elif CCALGO==STL_BST
/********************************************************************************
*** find coincidences - STL BST
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
        for (double htag : tmptags) {
            if (binary_search(ttags->begin(), ttags->end(), htag)) {
                cc_tags.emplace_back(htag);
            }
        }
        ccdata.cc_h = cc_tags.size();
        ccdata.cc_h_tags = cc_tags;
        
        // V
        tmptags = *vtags;
        std::ranges::for_each(tmptags.begin(), tmptags.end(), [&os](double &vt) { vt= round((ht-os)*10)/10; });
        cc_tags = {};
        for (double vtag : tmptags) {
            if (binary_search(ttags->begin(), ttags->end(), vtag)) {
                cc_tags.emplace_back(vtag);
            }
        }
        ccdata.cc_v = cc_tags.size();
        ccdata.cc_v_tags = cc_tags;
        ccdata.offset = os;
        result[i] = ccdata;
    }
    return result;
}
#endif

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
    if (TRUNCATE_S > 0) {
        string truncstr = string(3 - to_string(TRUNCATE_S).length(), '0') + to_string(TRUNCATE_S);
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
                HIST_START = stod(value);
            } else if (name == "hist_stop") {
                HIST_STOP = stod(value);
            } else if (name == "hist_step") {
                HIST_STEP = stod(value);
            } else if (name == "chn_h") {
                CHN_H = stoi(value);
            } else if (name == "chn_v") {
                CHN_V = stoi(value);
            } else if (name == "chn_tr") {
                CHN_TR = stoi(value);
            } else if (name == "chn_fpga") {
                CHN_FPGA = stoi(value);
            } else if (name == "fpga_delay_min") {
                FPGA_T_MIN = stod(value);
            } else if (name == "fpga_delay_max") {
                FPGA_T_MAX = stod(value);
            } else if (name == "truncate") {
                TRUNCATE_S = stoi(value);
            } else {
                cout << "unknown config name: " << name << endl;
            }
        }
        cfgf.close();
    } else {
        cerr << "Couldn't open config file for reading.\n";
    }
    
    cout << "---config---" << endl;
    cout << "hist_start : " << HIST_START << "ns" << endl;
    cout << "hist_stop  : " << HIST_STOP  << "ns" << endl;
    cout << "hist_step  : " << HIST_STEP  << "ns" << endl;
    cout << "thigger chn: " << CHN_TR     << endl;
    cout << "fpga chn   : " << CHN_FPGA   << endl;
    cout << "fpga t_min : " << FPGA_T_MIN << "ns" << endl;
    cout << "fpga t_max : " << FPGA_T_MAX << "ns" << endl;
    cout << "H chn      : " << CHN_H      << endl;
    cout << "V chn      : " << CHN_V      << endl;
    cout << "truncate at: " << TRUNCATE_S << "s" << endl;
}

