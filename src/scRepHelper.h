#include <Rcpp.h>
#include <vector>

class scRepHelper {
public: 
    static long double sum(std::vector<long double>& v) {
        long double n = 0;
        for (long double num : v) {
            n += num;
        }
        return n;
    }

    static Rcpp::NumericVector convertZerosToNA(std::vector<long double>& v, int len) {
        Rcpp::NumericVector converted (len, R_NaReal);
        for (int i = 0; i < len; i++) {
            if (v[i] > 0) {
                converted[i] = v[i];
            }
        }
        return converted;
    }
};