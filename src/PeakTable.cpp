//
// Created by Gregory Johnson on 8/22/25.
//

#include "PeakTable.h"

#include "Math.h"

PeakTable::PeakTable(const Rcpp::DataFrame& peakTable, const std::vector<std::string>& uniqueSampleList,
                     const double cvCutOff, const size_t replicates, bool fixPeaks) {
    const std::vector<std::string>& compounds = peakTable["Compound"];
    const std::vector<std::string>& sampleCodes = peakTable["Sample_Code"];
    const std::vector<std::string>& injections = peakTable["sample"];
    const std::vector<double> intensity = peakTable["intensity"];
    const size_t size = compounds.size();
    features = std::vector<FeatureData>(size);
    coefficientOfVariance = std::vector<std::list<double>>(replicates + 1);
    int index = 0;
    for (const auto& sample : uniqueSampleList) {
        sampleCodesToIndex[sample] = index++;
    }
    for (size_t i = 0, featureIndex = 0; i < size; i++) {
        const std::string& currentCompound = compounds[i];
        features[featureIndex].compoundName = currentCompound;
        features[featureIndex].metaData.intensityPerSample = std::vector<std::vector<double>>(uniqueSampleList.size());
        while (compounds[i] == currentCompound) {
            // fill in meta data, all the data should be sorted and we should
            // Get all the data for one compound as we iterate
            const int sampleIndex = sampleCodesToIndex[sampleCodes[i]];
            features[featureIndex].metaData.intensityPerSample[sampleIndex].emplace_back(intensity[i]);
            i++;
        }
        int idx = 0;
        for (const auto& intensityList : features[featureIndex].metaData.intensityPerSample) {
            bool hasPassedCv = false;
            compoundNamesToCV.emplace_back(currentCompound);
            sampleCodeList.emplace_back(uniqueSampleList[idx++]);
            double cvScore = VectorMath::CoefficientOfVarianceCalculation(intensityList);
            if (cvScore < cvCutOff && cvScore > 0) {
                hasPassedCv = true;
            }
            coefficientOfVariance[0].emplace_back(cvScore);
            if (fixPeaks) {
                for (size_t j = 0; j < replicates; j++) {
                    std::vector<double> permutations(replicates - 1);
                    size_t count = 0;
                    for (size_t k = 0; k < replicates; k++) {
                        if (j == k) continue;
                        permutations[count++] = intensityList[k];
                    }
                    cvScore = VectorMath::CoefficientOfVarianceCalculation(permutations);
                    if (cvScore < cvCutOff && cvScore > 0)
                        hasPassedCv = true;
                    coefficientOfVariance[j+1].emplace_back(cvScore);
                }
            }
            passesCV.emplace_back(hasPassedCv);
        }
        featureIndex++;
        i--;
    }

}

Rcpp::DataFrame PeakTable::GetCVTable() const {
    Rcpp::DataFrame df = Rcpp::DataFrame::create(
        Rcpp::Named("Compound") = compoundNamesToCV,
        Rcpp::Named("Sample_Code") = sampleCodeList,
        Rcpp::Named("PassesCvFilter") = passesCV);
    int i = 1;
    for (const auto& cvScores : coefficientOfVariance) {
        if (cvScores.size() <= 0) continue;
        df.push_back(cvScores, "CVPermutation_" + std::to_string(i++));
    }
    return df;
}
