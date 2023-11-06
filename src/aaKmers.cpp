// // naive implementation of hashmap-based AA counting
// // By Qile Yang

// #include <Rcpp.h>
// #include <vector>
// #include <string>
// #include "screpUtils.h"

// // // the following could be used if a 5-bit based AA kmer counting is needed
// // inline unsigned short int toAAIndex(const char aa) {
// //     switch(aa) {
// //         case 'A': return 0;
// //         case 'C': return 1;
// //         case 'D': return 2;
// //         case 'E': return 3;
// //         case 'F': return 4;
// //         case 'G': return 5;
// //         case 'H': return 6;
// //         case 'I': return 7;
// //         case 'K': return 8;
// //         case 'L': return 9;
// //         case 'M': return 10;
// //         case 'N': return 11;
// //         case 'P': return 12;
// //         case 'Q': return 13;
// //         case 'R': return 14;
// //         case 'S': return 15;
// //         case 'T': return 16;
// //         case 'V': return 17;
// //         case 'W': return 18;
// //         default: return 19; // case 'Y'
// //     }
// // }

// inline bool isAA(const char aa) {
//     switch(aa) {
//         case 'A': case 'C': case 'D': case 'E': case 'F': case 'G': case 'H': case 'I':
//         case 'K': case 'L': case 'M': case 'N': case 'P': case 'Q': case 'R': case 'S':
//         case 'T': case 'V': case 'W': case 'Y': return true;
//         default: return false;
//     } // may be better to input a set of aa's :P hopefully the compiler doesnt reconstruct the hashmap every time this function is called
// }

// inline void updateSkip(int& skip, const char c, const int k) {
//     if (!isAA(c)) {
//         skip = k;
//     } else if (skip > 0) {
//         skip--;
//     }
// }

// // [[Rcpp::export]]
// Rcpp::NumericVector rcppGetAAKmerPercent(
//     const std::vector<std::string>& seqs, const std::vector<std::string>& motifs, const int k
// ) {
//     std::unordered_map<std::string, int> map = toUnorderedMap(motifs);
//     std::vector<long double> counts (motifs.size(), 0);
//     int skip = k;

//     for (std::string seq : seqs) {
//         for (int i = 0; i < (int) seq.size(); i++) {
//             updateSkip(skip, seq[i], k);
//             if (skip > 0) {
//                 continue;
//             }
//             std::string kmer = seq.substr(i - k + 1, k);
//             if (map.find(kmer) != map.end()) {
//                 counts[map[kmer]]++;
//             }
//         }
//     }
//     return convertZerosToNA(counts, motifs.size());
// }