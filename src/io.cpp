#include "io.h"

/********************************************************************************
*** check if a file exists
*/
bool fileExists(const std::string& fn){
    struct stat buf;
    return (stat(fn.c_str(), &buf) != -1);
}

/********************************************************************************
*** read capnp tag file
*/
void readcapnptags(const std::string &fn, std::vector<long long> &result){
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = CAPNP_TRAVERSAL_LIMIT;
    if (std::filesystem::path(fn).extension() == ".zst") {
        std::stringstream ss;
        std::ifstream input(fn,std::ios_base::in);
        boost::iostreams::filtering_streambuf<boost::iostreams::input>bins;
        bins.push(boost::iostreams::zstd_decompressor());
        bins.push(input);
        boost::iostreams::copy(bins,ss);

        ::kj::std::StdInputStream stream(ss);
        try {
            while (true) {
                ::capnp::InputStreamMessageReader message(stream, opts);
                
                Tags::Reader reader = message.getRoot<Tags>();
                auto tagllist = reader.getTags();
                
                for (auto list: tagllist) {
                    for (auto tags: list) {
                        result.emplace_back(tags.getChannel());
                        result.emplace_back(tags.getTime());
                    }
                }
            }
        } catch (::kj::Exception &kje) {
            if (int(kje.getType())==0) {
                //expected EOF
            } else {
                std::cout << "kj exception of type " << int(kje.getType()) << std::endl;
                std::cout << std::string(kje.getDescription()) << std::endl;
            }
        } catch (std::exception &se) {
            std::cout << se.what() << std::endl;
        }
    } else {
        int fd = open(fn.c_str(), O_RDONLY);
        
        try {
            while (true) {
                ::capnp::StreamFdMessageReader message(fd, opts);
                Tags::Reader reader = message.getRoot<Tags>();
                auto tagllist = reader.getTags();
                for (auto list: tagllist) {
                    for (auto tags: list) {
                        result.emplace_back(tags.getChannel());
                        result.emplace_back(tags.getTime());
                    }
                }
            }
            close(fd);
        } catch (::kj::Exception &kje) {
            if (int(kje.getType())==0) {
                //expected EOF
            } else {
                std::cout << "kj exception of type " << int(kje.getType()) << std::endl;
                std::cout << std::string(kje.getDescription()) << std::endl;
            }
        } catch (std::exception &se) {
            std::cout << se.what() << std::endl;
        }
    }
}

void writecapnptags(std::string &fn, std::vector<long long> data, const bool compress, const uint8_t compression_level) {
    const size_t maxlistlen = (1<<28)-1;
    auto numlists = (uint32_t) (data.size()/maxlistlen);
    if (data.size()>maxlistlen*numlists) {
        ++numlists;
    }
    std::vector<std::vector<long long>> chunks;
    for (size_t i = 0; i<numlists; ++i) {
        auto end = ((i+1)*maxlistlen > data.size()) ? data.end() : data.begin()+(i+1)*maxlistlen;
        chunks.emplace_back(std::vector<long long>(data.begin()+i*maxlistlen, end));
    }
    
    if (compress) {
        fn = fn + ".zst";
        ::capnp::MallocMessageBuilder message;
        auto tags = message.initRoot<Tags>();
        auto llist = tags.initTags(numlists);
        for (size_t i=0; i<numlists; ++i){
            auto list = llist.init(i, chunks[i].size()/2);
            for (size_t j=0; j<chunks[i].size()/2; ++j){
                list[j].setChannel(chunks[i][2*j]);
                list[j].setTime(chunks[i][2*j+1]);
            }
        }
        std::ofstream ofs (fn, std::ios::out | std::ios::binary); 
        std::stringstream ss;
        ::kj::std::StdOutputStream kjos(ss);
        
        boost::iostreams::filtering_streambuf<boost::iostreams::output> outStream; 
        outStream.push(boost::iostreams::zstd_compressor(compression_level)); 
        outStream.push(ofs); 
        
        writeMessage(kjos, message);
        boost::iostreams::copy(ss, outStream); 
    } else {
        ::capnp::MallocMessageBuilder message;
        auto tags = message.initRoot<Tags>();
        auto llist = tags.initTags(numlists);
        for (size_t i=0; i<numlists; ++i){
            auto list = llist.init(i, chunks[i].size()/2);
            for (size_t j=0; j<chunks[i].size()/2; ++j){
                list[j].setChannel(chunks[i][2*j]);
                list[j].setTime(chunks[i][2*j+1]);
            }
        }

        int fd = open(fn.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
        writeMessageToFd(fd, message);
        close(fd);
    }

}

/********************************************************************************
*** read tsv tag file
*/
void readTSVtags(const std::string &fn, std::vector<long long> &result) {
    if (FILE *f = fopen(fn.c_str(), "r")) {
        (void)fseek(f, 0, SEEK_END);
        result.reserve(2*ftell(f));
        (void)fclose(f);
    }
    std::ifstream f(fn);
    std::stringstream ss;
    std::string s;
    ss << f.rdbuf();    
    f.close();
    long long int val = 0;
    while(std::getline(ss,s,'\t')) {
        std::from_chars(s.data(), s.data()+s.size(), val);
        result.push_back(val);
        std::getline(ss,s,'\n');
        std::from_chars(s.data(), s.data()+s.size(), val);
        result.push_back(val);
    }
}

void writeHDFtags(const std::string &fn, const std::vector<long long> &r, const uint8_t compression_alg, const uint8_t compression_level) {     
    hsize_t dim_values[2];
    dim_values[0] = size_t(r.size()/2);
    dim_values[1] = 2;
    int rank_values = sizeof(dim_values) / sizeof(hsize_t);
    //chunk for compression
    hsize_t chunkdims_values[2];
    if (r.size() > 16380) {
        chunkdims_values[0] = 16380;//what pandas uses
    } else {
        chunkdims_values[0] = uint(r.size()/2);
    }
    chunkdims_values[1] = 2;
    // group name to be compatible with python pandas
    std::string groupName = "/tags";
    std::string datasetPath_values = "/tags/block0_values";
    
    H5::H5File *file = new H5::H5File(fn.c_str(), H5F_ACC_TRUNC);
    H5::Group *group = new H5::Group(file->createGroup(groupName.c_str()));
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

void lltoTSV(const std::string &fn, const std::vector<long long> &data) {
    auto tsvfile = fmt::output_file(fn, fmt::buffer_size=262144);
    size_t len = data.size();
    for (size_t i = 0; i<len-1; i+=2) {
        tsvfile.print("{0}\t{1}\n", data[i], data[i+1]);
    }
}

void readHDF5tags(const std::string &fn, std::vector<long long>& result)
{
    std::string datasetPath = "/tags/block0_values";
    H5::H5File file(fn.c_str(), H5F_ACC_RDONLY);

    // access the dataset by path
    H5::DataSet dataset = file.openDataSet(datasetPath.c_str());

    // get the dataspace
    H5::DataSpace dataspace = dataset.getSpace();
  
    // Get the dimension of dataspace
    hsize_t dims_out[2];
    dataspace.getSimpleExtentDims( dims_out, nullptr);
    hsize_t NX=dims_out[0];
    hsize_t NY=dims_out[1];
    assert(NY==2);
    
    // read
    long long* data_out = new long long[NX*NY];
    
    long long out_data_len = NX*NY;
    dataset.read( data_out, H5::PredType::NATIVE_LLONG);
    result = std::vector<long long>(data_out, data_out+out_data_len);
    file.close();
    delete[] data_out;
}

std::vector<std::string> get_new_tagfiles(const std::string &raw_ext, const std::string &analyzed_ext, const std::vector<std::string> &excludes) {
    std::string path = std::filesystem::current_path();
    std::vector<std::string> alltagfiles = {};
    std::vector<std::string> analyzedfiles = {};
    std::vector<std::string> newtagfiles = {};
    
    for ( const auto & p : std::filesystem::directory_iterator(path) ) {
        if ( p.path().extension() == raw_ext ) {
            if (std::find(begin(excludes), end(excludes), p.path().filename()) == excludes.end()) { 
                alltagfiles.emplace_back(p.path());
            }
        } else if ((!analyzed_ext.empty()) && ( p.path().extension() == analyzed_ext )) {
            analyzedfiles.emplace_back(p.path());
        }
    }
    
    newtagfiles = alltagfiles;
    if ( analyzedfiles.empty() ) {
        return alltagfiles;
    }
    
    auto it = newtagfiles.begin();
    while (it != newtagfiles.end() ) {
        bool found = false;
        std::string analyzed_string = *it;
        stringreplace(analyzed_string, raw_ext, "");
        for ( const auto &p : analyzedfiles ) {
            std::string tmpstring = p;
            stringreplace(tmpstring,analyzed_ext,"");
            
            if (( tmpstring == analyzed_string)) {
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
