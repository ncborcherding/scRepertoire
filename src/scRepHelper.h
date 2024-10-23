#pragma once

#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>

#define BarcodeIndciesMap std::unordered_map<std::string, std::vector<int>>

class scRepHelper {
public:
    static double sum(std::vector<double>& v) {
        double n = 0;
        for (double num : v) {
            n += num;
        }
        return n;
    }

    static Rcpp::NumericVector convertZerosToNA(std::vector<double>& v, int len) {
        Rcpp::NumericVector converted (len, R_NaReal);
        for (int i = 0; i < len; i++) {
            if (v[i] > 0) {
                converted[i] = v[i];
            }
        }
        return converted;
    }

    static Rcpp::NumericVector convertZerosToNA(std::vector<double>& v) {
        return convertZerosToNA(v, (int) v.size());
    }

    static std::vector<std::vector<int>> constructBarcodeIndex(
        std::vector<std::string>& conDfBarcodes, std::vector<std::string> data2Barcodes
    ) {
        std::vector<std::vector<int>> outputBarcodeIndex (conDfBarcodes.size());
        BarcodeIndciesMap data2BarcodeIndiciesMap = stringToIndiciesMap(data2Barcodes);

        for (int i = 0; i < (int) conDfBarcodes.size(); i++) {
            std::string& barcode = conDfBarcodes[i];
            if (data2BarcodeIndiciesMap.find(barcode) != data2BarcodeIndiciesMap.end()) {
                outputBarcodeIndex[i] = data2BarcodeIndiciesMap[barcode];
            }
        }
        return outputBarcodeIndex;
    }

    static std::unordered_map<std::string, std::vector<int>> stringToIndiciesMap(
        std::vector<std::string>& v
    ) {
        std::unordered_map<std::string, std::vector<int>> map;
        for (int i = 0; i < (int) v.size(); i++) {
            map[v[i]].push_back(i);
        }
        return map;
    }

    static std::vector<std::vector<std::string>> initStringMatrix(
        int x, int y, std::string initVal
    ) {
        std::vector<std::vector<std::string>> stringMatrix (
            x, std::vector<std::string> (y, initVal)
        );
        return stringMatrix;
    }
};
