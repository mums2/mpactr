//
// Created by Gregory Johnson on 8/22/25.
//

#ifndef MPACTR_PEAKTABLE_H
#define MPACTR_PEAKTABLE_H
#include <vector>
#include <Rcpp.h>
#include "FeatureData.h"


class CvFilter {
public:
    explicit CvFilter(const bool useRecursion):useRecursiveMethod(useRecursion) {}
    void CalculateCV(const Rcpp::DataFrame& peakTable, const std::vector<std::string>& uniqueSampleList,
        double cvCutOff, size_t replicates);
    Rcpp::DataFrame GetCvTable() const;
private:
    std::vector<FeatureData> features;
    std::list<std::string> sampleCodeList;
    std::list<std::string> biologicalGroupsList;
    std::list<bool> passesCV;
    std::list<double> coefficientOfVariance;
    std::list<std::string> compoundNamesToCV;
    std::unordered_map<std::string, int> sampleCodesToIndex;
    bool useRecursiveMethod;
};


#endif //MPACTR_PEAKTABLE_H