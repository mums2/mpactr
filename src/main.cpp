#include <Rcpp.h>
#include <iostream>

// [[Rcpp::export]]
void Test() {
    Rcpp::Rcout << "Hello World\n";
}


void FilterMispickedIons(const Rcpp::DataFrame& peakTable, const double ringWin, const double isoWin,
    const double trWin, const double maxIsoShift, const bool mergePeaks,
    const std::string& mergeMethod) {
    const Rcpp::NumericVector mzVector = peakTable["mz"];
    const Rcpp::NumericVector rtVector = peakTable["rt"];
    const Rcpp::NumericVector compoundVector = peakTable["compound"];
    const auto rows = static_cast<size_t>(peakTable.nrows());

    for(size_t i = 0; i < rows - 1; i++) {
        
    }





}