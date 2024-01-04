// Rcpp replacement for .parseTCR
// By Qile Yang

#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_map>
#include "scRepHelper.h"

#define BarcodeIndciesMap std::unordered_map<std::string, std::vector<int>>

std::vector<std::vector<int>> constructBarcodeIndex(
    std::vector<std::string>& conDfBarcodes, std::vector<std::string> data2Barcodes
) {
    std::vector<std::vector<int>> outputBarcodeIndex (conDfBarcodes.size());
    BarcodeIndciesMap data2BarcodeIndiciesMap = scRepHelper::stringToIndiciesMap(data2Barcodes);

    for (int i = 0; i < (int) conDfBarcodes.size(); i++) {
        std::string& barcode = conDfBarcodes[i];
        if (data2BarcodeIndiciesMap.find(barcode) != data2BarcodeIndiciesMap.end()) {
            outputBarcodeIndex[i] = data2BarcodeIndiciesMap[barcode];
        }
    }
    return outputBarcodeIndex;
}

class TcrParser {
public:
    // variable for the eventual output Con.df
    std::vector<std::vector<std::string>> conDf;

    // variables for *references* to columns on data2
    Rcpp::CharacterVector data2ChainTypes;
    Rcpp::CharacterVector data2Tcr1;
    Rcpp::CharacterVector data2Tcr2;
    Rcpp::CharacterVector data2Cdr3;
    Rcpp::CharacterVector data2Cdr3Nt;

    // variable for helper barcode index
    std::vector<std::vector<int>> barcodeIndex;

    // constructor
    TcrParser(
        Rcpp::DataFrame& data2, std::vector<std::string>& uniqueData2Barcodes
    ) {
        // construct conDf, initializaing the matrix to "NA" *strings*
        conDf = scRepHelper::initStringMatrix(
            7, uniqueData2Barcodes.size(), "NA"
        );
        conDf[0] = uniqueData2Barcodes;

        // set references to fixed data2 columns
        data2ChainTypes = data2[data2.findName("chain")];
        data2Cdr3 = data2[data2.findName("cdr3")];
        data2Cdr3Nt = data2[data2.findName("cdr3_nt")];

        // setting reference to the TCR columns assuming all extra columns come before
        data2Tcr1 = data2[data2.findName("TCR1")];
        data2Tcr2 = data2[data2.findName("TCR2")];

        // construct barcodeIndex
        barcodeIndex = constructBarcodeIndex(
            uniqueData2Barcodes, Rcpp::as<std::vector<std::string>>(data2[data2.findName("barcode")])
        );
    }

    // Rcpp implementation of .parseTCR()
    void parseTCR() {
        for (int y = 0; y < (int) conDf[0].size(); y++) {
            for (int index : barcodeIndex[y]) {
                std::string chainType = std::string(data2ChainTypes[index]);
                if (chainType == "TRA" || chainType == "TRG") {
                    handleTcr1(y, index);
                } else if (chainType == "TRB" || chainType == "TRD") {
                    handleTcr2(y, index);
                }
            }
        }
    }

    // parseTCR() helpers

    void handleTcr1(int y, int data2index) {
        handleTcr(y, data2index, data2Tcr1, 1, 2, 3);
    }

    void handleTcr2(int y, int data2index) {
        handleTcr(y, data2index, data2Tcr2, 4, 5, 6);
    }

    void handleTcr(
        int y, int data2index, Rcpp::CharacterVector& data2tcr, int tcr, int cdr3aa, int cdr3nt
    ) {
        if (conDf[tcr][y] == "NA") { 
            conDf[tcr][y] = data2tcr[data2index];
            conDf[cdr3aa][y] = data2Cdr3[data2index];
            conDf[cdr3nt][y] = data2Cdr3Nt[data2index];
        } else {
            conDf[tcr][y] += ";" + data2tcr[data2index];
            conDf[cdr3aa][y] += ";" + data2Cdr3[data2index];
            conDf[cdr3nt][y] += ";" + data2Cdr3Nt[data2index];
        }
    }

    // return Con.df after TCR parsing
    Rcpp::DataFrame getConDf() {
        return Rcpp::DataFrame::create(
            Rcpp::Named("barcode") = conDf[0],
            Rcpp::Named("TCR1") = conDf[1],
            Rcpp::Named("cdr3_aa1") = conDf[2],
            Rcpp::Named("cdr3_nt1") = conDf[3],
            Rcpp::Named("TCR2") = conDf[4],
            Rcpp::Named("cdr3_aa2") = conDf[5],
            Rcpp::Named("cdr3_nt2") = conDf[6]
        );
    }
};

// [[Rcpp::export]]
Rcpp::DataFrame rcppConstructConDfAndParseTCR(
    Rcpp::DataFrame& data2, std::vector<std::string> uniqueData2Barcodes
) {
    TcrParser parser = TcrParser(data2, uniqueData2Barcodes);
    parser.parseTCR();
    return parser.getConDf();
}
