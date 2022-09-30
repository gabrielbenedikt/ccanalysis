#include "io_protobuf.h"

int ser_histogram_protobuf(const histograms* data, const std::string &fname) {
    std::cout << "writing to file " << fname << std::endl;
    histogramset::histograms hdat;
    hdat.set_resolution(data->resolution);
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
    }
    
    std::fstream ofs(fname, std::ios::out | std::ios::trunc | std::ios::binary);
    if (!hdat.SerializeToOstream(&ofs)) {
        std::cerr << "Failed to write." << std::endl;
        return -1;
    }
    ofs.close();
    
    return 0;
}

histograms deser_histogram_protobuf(const std::string &fn) {
    histograms h;
    histogramset::histograms pbhist;
    std::fstream ifs(fn, std::ios::in | std::ios::binary);

    if (!pbhist.ParseFromIstream(&ifs)) {
        std::cerr << "Failed to parse protobuf file." << std::endl;
    } else {
        h.resolution = pbhist.resolution();
        
        for (int i=0; i<pbhist.meastime_size(); ++i) {
            h.meastime.emplace_back(pbhist.meastime(i));
            
            std::vector<long long>              offsets;
            std::vector<uint32_t>               cc;
            std::vector<uint16_t>               pattern;
            std::vector<std::vector<long long>> cc_tags;
            
            const auto &pb_offsets      = pbhist.offsets(i);
            const auto &pb_cc           = pbhist.cc(i);
            const auto &pb_pattern      = pbhist.pattern(i);
            const auto &pb_cc_tags_v    = pbhist.cc_tags(i);
            
            offsets.reserve(pb_offsets.arr_size());
            cc.reserve(pb_cc.arr_size());
            pattern.reserve(pb_pattern.arr_size());
            cc_tags.reserve(pb_cc_tags_v.arr_size());
            
            for (const auto &pbo: pb_offsets.arr()) {
               offsets.emplace_back(pbo);
            }
            for (const auto &pbcc: pb_cc.arr()) {
               cc.emplace_back(pbcc);
            }
            for (const auto &pbp: pb_pattern.arr()) {
               pattern.emplace_back(pbp);
            }
            for (const auto &pbctv: pb_cc_tags_v.arr()) {
                std::vector<long long> ct;
                ct.reserve(pbctv.arr_size());
                for (const auto &pbct: pbctv.arr()) {
                    ct.emplace_back(pbct);
                }
                cc_tags.emplace_back(ct);
            }
            
            h.offsets.emplace_back(offsets);
            h.cc.emplace_back(cc);
            h.pattern.emplace_back(pattern);
            h.cc_tags.emplace_back(cc_tags);
        }
    }
    
    return h;
}
