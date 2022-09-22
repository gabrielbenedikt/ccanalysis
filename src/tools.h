#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <regex>
#include <iomanip>

bool stringreplace(std::string& str, const std::string& from, const std::string& with);
std::vector<std::vector<uint16_t>> parse_patterns(const std::string &instring);

/********************************************************************************
*** create string of integer, with fixed width
*/
template<typename T, typename U>
requires std::integral<T> && std::unsigned_integral<U>
std::string fixed_width_intstr(const T number, const U width){
    std::ostringstream ss;
    if (number < 0) {
        ss << '-';
    }
    ss << std::setfill('0') << std::setw(width) << (number < 0 ? -number : number);
    return ss.str();
}

/********************************************************************************
*** print elements of vector
*/
template<typename T>
void print_vector(const std::vector<T> v) {
    std::cout << "[";
    for (T e : v) {
        if (e != v.back()) {
            std::cout << e << ", ";
        } else {
            std::cout << e << "]" << std::endl;
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
        if (std::find(result.begin(), result.end(), e) == result.end()) {
            result.push_back(e);
        }
    }
    return result;
}

/********************************************************************************
*** create vector holding range of values
*/
template<class T, typename std::enable_if<std::is_floating_point<T>::value>::type...>
std::vector<T> arange(const T start, const T stop, const T step) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step) {
        values.emplace_back(value);
    }
    return values;
}

template<class T, typename std::enable_if<std::is_integral<T>::value>::type...>
std::vector<T> arange(const T start, const T stop, const T step) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step) {
        values.emplace_back(value);
    }
    return values;
}

template<class T>
requires std::integral<T>
T roundto(const T num, const T to)
{
    if (to == 0) {
        return num;
    }
    
    T r = abs(num) % to;
    if (r == 0) {
        return num;
    }
    
    if (num < 0) {
        if (2*r > to) {
            return num + r - to;
        } else {
            return num + r;
        }
    } else {
        if (2*r > to) {
            return num - r + to;
        } else {
            return num - r;
        }
    }
}
