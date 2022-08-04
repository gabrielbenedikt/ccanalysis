#include "io.h"

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
*** read tsv tag file
*/
void readTSVtags(const std::string fn, std::vector<long long> &result, long long& out_data_len) {
    if (FILE *f = fopen(fn.c_str(), "r")) {
        
        fseek(f, 0, SEEK_END);
        long long fsize = ftell(f);
        fseek(f, 0, SEEK_SET);
        result.reserve(fsize);
        
        const uint16_t MAX_LINESIZE = 256;
        char line[MAX_LINESIZE];
        const char *delim = "\t";
        while(fgets(line, MAX_LINESIZE, f) != NULL) {
            result.push_back(std::stoll(strtok(line, delim)));
            result.push_back(std::stoll(strtok(NULL, delim)));;
        }
    }
    
    out_data_len = result.size();
}

void writeHDFtags(const std::string fn, const std::vector<long long> &r, const uint8_t compression_alg, const uint8_t compression_level) {     
    hsize_t dim_values[2];
    dim_values[0] = size_t(r.size()/2);
    dim_values[1] = 2;
    int rank_values = sizeof(dim_values) / sizeof(hsize_t);
    //chunk for compression
    hsize_t chunkdims_values[2];
    chunkdims_values[0] = 16380;//what pandas uses
    chunkdims_values[1] = 2;
    // group name to be compatible with python pandas
    std::string groupName = "/tags";
    std::string datasetPath_values = "/tags/block0_values";
    
    H5::H5File *file = new H5::H5File(fn.c_str(), H5F_ACC_TRUNC);
    H5::Group* group = new H5::Group(file->createGroup(groupName.c_str()));
    H5::DSetCreatPropList ds_creat_plist_values;
    ds_creat_plist_values.setChunk(2, chunkdims_values);
    
    switch (compression_alg) {
        case HDF5_COMP_ALG_ZLIB:
            ds_creat_plist_values.setShuffle();
            ds_creat_plist_values.setDeflate(compression_level);
            break;
        case HDF5_COMP_ALG_SZIP:
            ds_creat_plist_values.setShuffle();
            ds_creat_plist_values.setSzip(H5_SZIP_NN_OPTION_MASK, compression_level);
            break;
        case HDF5_COMP_ALG_NBIT:
            ds_creat_plist_values.setShuffle();
            ds_creat_plist_values.setNbit();
            break;
        default:
            break;
    }
    ds_creat_plist_values.setFillTime(H5D_FILL_TIME_ALLOC);
    ds_creat_plist_values.setAllocTime(H5D_ALLOC_TIME_INCR);
    
    H5::DataSpace space_values(rank_values, dim_values);
    H5::DataSet* dataset_values = new H5::DataSet(file->createDataSet(datasetPath_values.c_str(), H5::PredType::NATIVE_LLONG, space_values, ds_creat_plist_values));
    dataset_values->write(r.data(), H5::PredType::NATIVE_LLONG);
    dataset_values->close();
    delete dataset_values;
    delete group;
    delete file;
}

void writeHDFtagsC(const std::string fn, const std::vector<long long> &r, const uint8_t compression_alg, const uint8_t compression_level) { 
    hid_t file_id;
    hid_t ds_create_plist_values;
    hid_t dataset_id_values;
    hid_t dataspace_id_values;
    hid_t group_id;
    
    std::string groupName = "/tags";
    std::string datasetPath_values = "/tags/block0_values";
    
    //define data dimensions
    hsize_t dim_values[2];
    dim_values[0] = size_t(r.size()/2);
    dim_values[1] = 2;
    int rank_values = sizeof(dim_values) / sizeof(hsize_t);
    //chunk for compression
    hsize_t chunkdims_values[2];
    chunkdims_values[0] = 16380;//what pandas uses
    chunkdims_values[1] = 2;
    
    //oen file
    file_id = H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    // define dataset properties
    // group name to be compatible with python pandas
    ds_create_plist_values = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(ds_create_plist_values, 2, chunkdims_values);
    switch (compression_alg) {
        case HDF5_COMP_ALG_ZLIB:
            H5Pset_shuffle(ds_create_plist_values);
            H5Pset_deflate(ds_create_plist_values, compression_level);
            break;
        case HDF5_COMP_ALG_SZIP:
            H5Pset_shuffle(ds_create_plist_values);
            H5Pset_deflate(H5_SZIP_NN_OPTION_MASK, compression_level);
            break;
        case HDF5_COMP_ALG_NBIT:
            H5Pset_shuffle(ds_create_plist_values);
            H5Pset_nbit(ds_create_plist_values);
            break;
        default:
            break;
    }
    
    H5Pset_fill_time(ds_create_plist_values, H5D_FILL_TIME_ALLOC);
    H5Pset_alloc_time(ds_create_plist_values, H5D_ALLOC_TIME_INCR);
    
    //dataspaces
    dataspace_id_values = H5Screate_simple(rank_values, dim_values, NULL);
    
    //create group
    group_id = H5Gcreate2(file_id, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // datasets
    dataset_id_values   = H5Dcreate2(file_id, datasetPath_values.c_str(),   H5T_STD_I64LE,  dataspace_id_values,    H5P_DEFAULT, ds_create_plist_values, H5P_DEFAULT);
    
    //
    // write data
    //
    H5Dwrite(dataset_id_values, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, r.data());
    // close everything
    H5Dclose(dataset_id_values);
    H5Sclose(dataspace_id_values);
    H5Gclose(group_id);
    H5Fclose(file_id);
}

void lltoTSV(const std::string fn, const long long* data, const long long len) {
    std::ofstream tsvfile(fn);
    for (long long i = 0; i<len-1; i+=2) {
        tsvfile << data[i] << "\t" << data[i+1] << "\n";
    }
}

long long* readHDF5tags(const std::string fn, long long& out_data_len)
{
    std::string datasetPath = "/tags/block0_values";
    H5::H5File file(fn.c_str(), H5F_ACC_RDONLY);

    // access the dataset by path
    H5::DataSet dataset = file.openDataSet(datasetPath.c_str());

    // get the dataspace
    H5::DataSpace dataspace = dataset.getSpace();
  
    // Get the dimension of dataspace
    hsize_t dims_out[2];
    dataspace.getSimpleExtentDims( dims_out, NULL);
    int NX=dims_out[0];
    int NY=dims_out[1];
    assert(NY==2);
    
    // read
    long long* data_out = new long long[NX*NY];
    out_data_len = NX*NY;
    dataset.read( data_out, H5::PredType::NATIVE_LLONG);
    
    file.close();
    
    return data_out;
}


std::vector<std::string> get_new_tagfiles(const std::string raw_ext, const std::string analyzed_ext, const std::vector<std::string> excludes) {
    std::string path = std::filesystem::current_path();
    std::vector<std::string> alltagfiles = {};
    std::vector<std::string> analyzedfiles = {};
    std::vector<std::string> newtagfiles = {};
    
    for ( auto & p : std::filesystem::directory_iterator(path) ) {
        std::cout << "found file " << p << std::endl;
        if ( p.path().extension() == raw_ext ) {
            std::cout << "raw_ext found" << std::endl;
            if (std::find(begin(excludes), end(excludes), p.path().filename()) == excludes.end()) { 
                std::cout << "insert in alltagfiles" << std::endl;
                alltagfiles.emplace_back(p.path());
            }
        } else if ((analyzed_ext != "") && ( p.path().extension() == analyzed_ext )) {
            std::cout << "analyzed_ext found. insert in analyzed files" << std::endl;
            analyzedfiles.emplace_back(p.path());
        }
    }
    
    newtagfiles = alltagfiles;
    if ( analyzedfiles.size() == 0 ) {
        return alltagfiles;
    } else {
        auto it = newtagfiles.begin();
        while (it != newtagfiles.end() ) {
            bool found = false;
            std::string analyzed_string = *it;
            stringreplace(analyzed_string, raw_ext, "");
            for ( auto p : analyzedfiles ) {
                std::string tmpstring = p;
                stringreplace(tmpstring,analyzed_ext,"");
                
                if (( tmpstring == analyzed_string)) {
                    it=newtagfiles.erase(it);
                    found = true;
                    std::cout << p << " in analyzed files. remove for list to process" << std::endl;
                }
            }
            if (!found) {
                ++it;
            }
        }
        
        return newtagfiles;
    }
}