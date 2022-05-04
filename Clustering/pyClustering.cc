#include <vector>
#include "ClusteringAlgorithm.hh"

#include <pybind11.h>
#include <stl.h>

// #include <pybind11/numpy.h>÷
namespace py = pybind11;

using namespace std;

std::vector<int> clustering(int numEvt, std::vector<int> eventsListDirect, std::vector<int> copyListDirect, std::vector<int> strandListDirect, std::vector<int> eventsListIndirect, std::vector<int> copyListIndirect, std::vector<int> strandListIndirect);

std::vector<int> clustering(int numEvt, std::vector<int> eventsListDirect, std::vector<int> copyListDirect, std::vector<int> strandListDirect, std::vector<int> eventsListIndirect, std::vector<int> copyListIndirect, std::vector<int> strandListIndirect)
{
    std::vector<int> results(12, 0);

    // start clustering object
    ClusteringAlgorithm *fpClusteringDirect = new ClusteringAlgorithm(10, 2, 5, 37.5);   // eV
    ClusteringAlgorithm *fpClusteringIndirect = new ClusteringAlgorithm(10, 2, 5, 37.5); // eV
    ClusteringAlgorithm *fpClusteringTotal = new ClusteringAlgorithm(10, 2, 5, 37.5);    // eV

    // Clustering
    for (int e = 0; e < numEvt; ++e)
    {
        for (int i = 0; i < eventsListDirect.size(); ++i)
        {
            if (e == eventsListDirect[i])
            {
                fpClusteringDirect->RegisterDamage(copyListDirect[i], 0, strandListDirect[i]);
                fpClusteringTotal->RegisterDamage(copyListDirect[i], 0, strandListDirect[i]); // edep not used
            }
        }

        for (int i = 0; i < eventsListIndirect.size(); ++i)
        {
            if (e == eventsListIndirect[i])
            {
                fpClusteringIndirect->RegisterDamage(copyListIndirect[i], 0, strandListIndirect[i]);
                fpClusteringTotal->RegisterDamage(copyListIndirect[i], 0, strandListIndirect[i]); // edep not used
            }
        }

        std::map<int, int> sizeDistributionDirect = fpClusteringDirect->RunClustering();

        results[0] += fpClusteringDirect->GetTotalSB();
        results[1] += fpClusteringDirect->GetSSB();
        results[2] += fpClusteringDirect->GetComplexSSB();
        results[3] += fpClusteringDirect->GetDSB();

        fpClusteringDirect->Purge();
        sizeDistributionDirect.clear();

        std::map<int, int> sizeDistributionIndirect = fpClusteringIndirect->RunClustering();

        results[4] += fpClusteringIndirect->GetTotalSB();
        results[5] += fpClusteringIndirect->GetSSB();
        results[6] += fpClusteringIndirect->GetComplexSSB();
        results[7] += fpClusteringIndirect->GetDSB();

        fpClusteringIndirect->Purge();
        sizeDistributionIndirect.clear();

        std::map<int, int> sizeDistribution = fpClusteringTotal->RunClustering();

        results[8] += fpClusteringTotal->GetTotalSB();
        results[9] += fpClusteringTotal->GetSSB();
        results[10] += fpClusteringTotal->GetComplexSSB();
        results[11] += fpClusteringTotal->GetDSB();

        fpClusteringTotal->Purge();

        sizeDistribution.clear();
    }

    delete fpClusteringDirect;
    delete fpClusteringIndirect;
    delete fpClusteringTotal;

    return results;
}

PYBIND11_MODULE(clustering, m)
{
    m.doc() = "clustering"; // optional module docstring

    m.def("clustering", &clustering, "A function that performs clustering");
}