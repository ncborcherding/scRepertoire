// naive implementation of hashmap-based AA counting, with a 5-bit AA representation.
// By Qile Yang

#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>
#include "scRepHelper.h"

std::unordered_map<char, unsigned long int> allAaMap() {
    std::unordered_map<char, unsigned long int> map;
    map['A'] = 0;
    map['C'] = 1;
    map['D'] = 2;
    map['E'] = 3;
    map['F'] = 4;
    map['G'] = 5;
    map['H'] = 6;
    map['I'] = 7;
    map['K'] = 8;
    map['L'] = 9;
    map['M'] = 10;
    map['N'] = 11;
    map['P'] = 12;
    map['Q'] = 13;
    map['R'] = 14;
    map['S'] = 15;
    map['T'] = 16;
    map['V'] = 17;
    map['W'] = 18;
    map['Y'] = 19;
    return map;
}

class AaKmerCounter {
public:
    // ideally these are all constants except bins
    std::unordered_map<unsigned long int, int> aaUIntKmerMap;
    int k;
    unsigned long int mask;
    std::unordered_map<char, unsigned long int> aaIndexMap;

    std::vector<long double> bins;

    // constructor
    AaKmerCounter(const std::vector<std::string>& motifs, const int _k) {
        aaIndexMap = allAaMap();
        k = _k;
        mask = (unsigned long int) ((1 << (_k * 5)) - 1);
        aaUIntKmerMap = toAaUIntKmerMap(motifs);
        bins = std::vector<long double> (motifs.size(), 0.0);
    }

    std::unordered_map<unsigned long int, int> toAaUIntKmerMap(const std::vector<std::string>& motifs) {
        std::unordered_map<unsigned long int, int> map;
        for (int i = 0; i < (int) motifs.size(); i++) {
            unsigned long int kmer = 0;
            for (char aa : motifs[i]) {
                kmer = (kmer << 5) | toAaIndex(aa);
            }
            map[kmer] = i;
        }
        return map;
    }

    inline unsigned long int toAaIndex(const char aa) {
        if (aaIndexMap.find(aa) == aaIndexMap.end()) {
            return 20;
        }
        return aaIndexMap[aa];
    }

    inline void updateSkip(int& skip, const char c) {
        if (aaIndexMap.find(c) == aaIndexMap.end()) {
            skip = k;
        } else if (skip > 0) {
            skip--;
        }  
    }

    void countKmers(const std::vector<std::string>& seqs) {
        for (std::string seq : seqs) {
            if ((int) seq.size() < k) {
                continue;
            }

            int skip = 0;
            unsigned long int kmer = 0;

            for (int i = 0; i < (k - 1); i++) { // this segment to initialize the kmer should be deletable if skip is adjusted?
                kmer = (kmer << 5) | toAaIndex(seq[i]);
                updateSkip(skip, seq[i]);
            }

            for (int i = (k - 1); i < (int) seq.size(); i++) {
                kmer = ((kmer << 5) & mask) | toAaIndex(seq[i]);
                updateSkip(skip, seq[i]);
                if (skip == 0) {
                    bins[aaUIntKmerMap[kmer]]++;
                }
            }
        }
    }

    std::vector<long double> getCounts() {
        return bins;
    }
};

// [[Rcpp::export]]
Rcpp::NumericVector rcppGetAaKmerPercent(
    const std::vector<std::string>& seqs, const std::vector<std::string>& motifs, const int k
) {

    AaKmerCounter counter = AaKmerCounter(motifs, k);
    counter.countKmers(seqs);
    std::vector<long double> bins = counter.getCounts();

    long double binSum = scRepHelper::sum(bins);
    if (binSum == 0.0) { // pretty sure this can only happen if there arent valid seqs?
        return Rcpp::NumericVector (motifs.size(), R_NaReal);
    }
    
    double scaleFactor = 1 / binSum;
    for (int i = 0; i < (int) motifs.size(); i++) {
        bins[i] *= scaleFactor;
    }

    return scRepHelper::convertZerosToNA(bins, motifs.size());
}
