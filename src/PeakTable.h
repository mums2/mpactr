//
// Created by Gregory Johnson on 8/22/25.
//

#ifndef MPACTR_PEAKTABLE_H
#define MPACTR_PEAKTABLE_H
#include <vector>
#include <Rcpp.h>
#include "FeatureData.h"


class PeakTable {
public:
    PeakTable(const Rcpp::DataFrame& peakTable, const std::vector<std::string>& uniqueSampleList,
        size_t replicates, bool fixPeaks);
    Rcpp::DataFrame GetCVTable() const;
private:
    std::vector<FeatureData> features;
    std::vector<std::list<double>> coefficientOfVariance;
    std::list<std::string> compoundNamesToCV;
    std::unordered_map<std::string, int> sampleCodesToIndex;
};


#endif //MPACTR_PEAKTABLE_H