#ifndef TOOLS_H
#define TOOLS_H

#include <string>
#include <vector>
#include <iostream>
#include <regex>

bool stringreplace(std::string& str, const std::string& from, const std::string& to);
std::vector<std::vector<uint16_t>> parse_patterns(const std::string instring);

/********************************************************************************
*** print elements of vector
*/
template<typename T>
void print_vector(const std::vector<T> v) {
    for (T e : v) {
        if (e != v.back()) {
            std::cout << e << ", ";
        } else {
            std::cout << e << std::endl;
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

/********************************************************************************
*** create vector holding range of values
*/
template<class T, std::enable_if<std::is_floating_point<T>::value>::type...>
std::vector<T> arange(const T start, const T stop, const T step) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.emplace_back(value);
    return values;
}

template<class T, std::enable_if<std::is_integral<T>::value>::type...>
std::vector<T> arange(const T start, const T stop, const T step) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.emplace_back(value);
    return values;
}

#endif //TOOLS_H
