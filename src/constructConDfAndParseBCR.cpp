// Upcoming Rcpp replacement for .parseBCR
// By Qile Yang

#include <Rcpp.h>
#include <string>
#include <vector>
#include "scRepHelper.h"

#define BarcodeIndciesMap std::unordered_map<std::string, std::vector<int>>

// Temporary function to improve theoretical time complexity of the R function
// remove when the script is done
// [[Rcpp::export]]
std::vector<std::vector<int>> rcppConstructBarcodeIndex(
    std::vector<std::string>& conDfBarcodes, std::vector<std::string> data2Barcodes
) {
    std::vector<std::vector<int>> outputBarcodeIndex (conDfBarcodes.size());
    BarcodeIndciesMap data2BarcodeIndiciesMap = scRepHelper::stringToRIndiciesMap(data2Barcodes);

    for (int i = 0; i < (int) conDfBarcodes.size(); i++) {
        std::string& barcode = conDfBarcodes[i];
        if (data2BarcodeIndiciesMap.find(barcode) != data2BarcodeIndiciesMap.end()) {
            outputBarcodeIndex[i] = data2BarcodeIndiciesMap[barcode];
        }
    }
    return outputBarcodeIndex;
}

// class BcrParser {
// public:

// };
