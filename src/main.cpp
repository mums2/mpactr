#include <Rcpp.h>
#include <iostream>
#include <vector>

// [[Rcpp::export]]
Rcpp::List FilterMispickedIons(const Rcpp::DataFrame& peakTable, const double ringWin, const double isoWin,
    const double trWin, const double maxIsoShift) {
    const auto rows = static_cast<size_t>(peakTable.nrows());
    const Rcpp::NumericVector mzVector = peakTable["mz"];
    const Rcpp::NumericVector rtVector = peakTable["rt"];
    const Rcpp::CharacterVector compoundVector = peakTable["Compound"];
    std::vector<std::string> cutIons;
    cutIons.reserve(rows);
    std::vector<bool> cutIonsChecked(rows, false);
    std::unordered_map<std::string, size_t> cutIonDict;
    std::vector<std::vector<std::string>> mergeGroupList(rows);
    std::vector<std::string> nameVector;
    std::vector<size_t> indexSwap(rows);
    nameVector.reserve(rows);
    mergeGroupList.reserve(rows);
    for(size_t i = 0; i < rows - 1; i++) {
        const double mz = mzVector[i];
        const double rt = rtVector[i];
        const Rcpp::String compound = compoundVector[i];
        bool addName = false;
        for (size_t j = i + 1; j < rows; j++) {
            const Rcpp::String nextCompound = compoundVector[j];
            if (cutIonsChecked[j]) continue;
            const double nextMz = mzVector[j];
            const double nextRt = rtVector[j];
            const double massDifference = nextMz - mz;
            const double kmdDifference = massDifference - std::floor(massDifference);
            const bool shiftDifference = std::abs(massDifference) > maxIsoShift - 0.4;
            const double rtDifference = nextRt - rt;
            const double ringBand = std::fmod(std::floor(std::abs(massDifference) * (1/ringWin)), (1/ringWin));
            const double doubleBand = kmdDifference - 0.5004 - (std::floor(massDifference) * 0.007);
            const bool betweenIonCalculation = std::abs(rtDifference) <= trWin && (massDifference <= maxIsoShift - 0.4) &&
                (ringBand ==0 || doubleBand < isoWin);

            if (!shiftDifference && betweenIonCalculation) {
                cutIons.emplace_back(nextCompound);
                cutIonsChecked[j] = true;
                mergeGroupList[i].emplace_back(nextCompound);
                addName = true;
            }
        }
        if (addName)
            nameVector.emplace_back(compound);
    }
    size_t size = nameVector.size();
    Rcpp::List mergeGroups(static_cast<int>(size));
    size_t count = 0;
    for (size_t i = 0; i < rows; i++) {
        if (mergeGroupList[i].empty()) continue;
        mergeGroups[count++] = mergeGroupList[i];
    }
    mergeGroups.attr("names") = nameVector;


    return Rcpp::List::create(Rcpp::Named("cut_ions") = cutIons,
        Rcpp::Named("merge_groups") = mergeGroups);
}
