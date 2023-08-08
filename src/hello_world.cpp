#include <Rcpp.h>

// [[Rcpp::export]]
void rcpp_hello_world() {
	Rcpp::Rcout << "hello, world!" << "\n";
}
