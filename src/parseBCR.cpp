// Upcoming Rcpp replacement for .parseBCR
// By Qile Yang

#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_map>
#include "scRepHelper.h"

Rcpp::DataFrame rcppParseBCR(Rcpp::DataFrame conDf, st uniqueDf, Rcpp::DataFrame data2) {
    BarcodeIndciesMap barcodeIndex = scRepHelper::constructBarcodeIndex(
        Rcpp::as<std::vector<std::string>>(uniqueDf.findName("")),
        Rcpp::as<std::vector<std::string>>(data2[data2.findName("barcode")])
    );
}
