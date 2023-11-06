#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>

const long double sum(std::vector<long double>& v) {
    long double n = 0;
    for (long double num : v) {
        n += num;
    }
    return n;
}

Rcpp::NumericVector convertZerosToNA(std::vector<long double>& v, int len) {
    Rcpp::NumericVector converted (len, R_NaReal);
    for (int i = 0; i < len; i++) {
        if (v[i] > 0) {
            converted[i] = v[i];
        }
    }
    return converted;
}

std::unordered_map<std::string, int> toUnorderedMap(std::vector<std::string>& motifs) {
    std::unordered_map<std::string, int> map;
    for (int i = 0; i < (int) motifs.size(); i++) {
        map[motifs[i]] = i;
    }
    return map;
}