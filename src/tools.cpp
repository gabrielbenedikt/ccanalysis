#include "tools.h"
/********************************************************************************
*** find and replace in a string (find 1, NOT N)
*/
bool stringreplace(std::string& str, const std::string& from, const std::string& with) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos) {
        return false;
    }
    str.replace(start_pos, from.length(), with);
    return true;
}

/********************************************************************************
*** parse string containing patternspatterns
*/
std::vector<std::vector<uint16_t>> parse_patterns(const std::string &instring) {
    std::string patstr = instring;
    if (!(patstr.starts_with("{{"))) {
        std::cout << "pattern string malformatted: " << patstr << std::endl;
    }
    if (!(patstr.ends_with("}}"))) {
        std::cout << "pattern string malformatted: " << patstr << std::endl;
    }
    
    patstr.erase(0,1);
    patstr.erase(patstr.length()-1,1);
    
    size_t vec_start_idx = 0;
    size_t vec_stop_idx = 0;
    std::vector<size_t> sep_idxs = {};
    size_t curr_sep_idx = 0;
    
    std::string vector_string;
    std::vector<uint16_t> pattern = {};
    std::vector<std::vector<uint16_t>> patterns = {};
    std::string channelstr;
    // separate patterns
    while (vec_stop_idx < patstr.length()-1) {
        vec_start_idx =  patstr.find('{', vec_start_idx);
        vec_stop_idx = patstr.find('}', vec_start_idx);
        vector_string = patstr.substr(vec_start_idx, vec_stop_idx-vec_start_idx+1);
        
        sep_idxs = {};
        sep_idxs.push_back(vector_string.find('{'));
        curr_sep_idx = vector_string.find(',', 0);
        while(curr_sep_idx != std::string::npos) {
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

