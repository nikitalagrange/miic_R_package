#ifndef MIIC_COMPUTE_INFO_H_
#define MIIC_COMPUTE_INFO_H_

#include "computation_cache.h"
#include "structure.h"

namespace miic {
namespace computation {

double* getAllInfoNEW(int* ptrAllData, const std::vector<int>& ptrAllLevels,
    const std::vector<int>& ptrVarIdx, int nbrUi, int* ptrZiIdx, int nbrZi,
    int ziPos, int sampleSize, int sampleSizeEff, int modCplx, int k23,
    double* looklog, structure::MemorySpace* memory,
    const std::vector<double>& weights, double** freqs1, bool test_mar,
    std::shared_ptr<CtermCache> cache);
double lnfactorial(int n, double* looklog);
double logchoose(int n, int k, double* looklog);
double logchoose(int n, int k, double* looklog, double** lookchoose);

}  // namespace computation
}  // namespace miic

#endif  // MIIC_COMPUTE_INFO_H_
