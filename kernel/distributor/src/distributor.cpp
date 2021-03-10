#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <vector>

#include <libmuscle/libmuscle.hpp>
#include <ymmsl/ymmsl.hpp>

#include <util/logger.h>
#include <util/typeNames.hpp>

using libmuscle::Data;
using libmuscle::DataConstRef;
using libmuscle::Instance;
using libmuscle::Message;
using ymmsl::Operator;

using util::logger;

class Distributor {
private:
    double t_cur;

    void computeAndSendIncrementalBoundaries(const std::vector<int32_t>& domain, std::vector<int32_t>& previousSolid);

    Instance& instance;
    const std::string outPath;

    std::vector<int32_t> previousSolid;
    // std::vector<std::vector<std::vector<unsigned char> > > aggregateDomain;
    std::vector<int32_t> aggregateDomain;
    std::vector<size_t> boxSize;

public:

    void execute();

    Distributor (Instance& instance_):
    instance(instance_),
    outPath (instance.get_setting_as<std::string> ("muscle_data_out_dir")) {
        const std::string logName = instance.get_setting_as<std::string> ("log_file");
        if(!logName.empty())
            logger().setFile(outPath + logName);

        previousSolid.clear();
        aggregateDomain.clear();
        boxSize = {0, 0, 0};
        logger() << "Initializing Distributor" << std::endl;
    }

    bool isSurfaceNode(int i, int j, int k, const std::vector<int32_t>& nodeTypes) const {
        return (nodeTypes[i*boxSize[Y]*boxSize[Z] + j*boxSize[Z] + k] != FLUID && hasIndirectNeighborWithValue(nodeTypes, i, j, k, FLUID));
    }

    bool hasIndirectNeighborWithValue(const std::vector<int32_t>& domain, int i, int j, int k, int32_t value) const {
        for (int ii = i - 1; ii <= i + 1; ii++) {
            if (ii == -1 || ii == boxSize[X]) continue;
            for (int jj = j - 1; jj <= j + 1; jj++) {
                if (jj == -1 || jj == boxSize[Y]) continue;
                for (int kk = k - 1; kk <= k + 1; kk++) {
                    if (kk == -1 || kk == boxSize[Z] || (ii == i && jj == j && kk == k)) continue;
                    if (domain[ii*boxSize[Y]*boxSize[Z] + jj*boxSize[Z] + kk] == value) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
};

void Distributor::execute() {
    logger() << "Receiving domain..." << std::endl;

    auto domainObs = instance.receive("domainIn");
    t_cur = domainObs.timestamp();
    logger() << "Domain received for time " << t_cur << std::endl;

    size_t doSize = domainObs.data().size();
    boxSize = domainObs.data().shape();
    logger() << "Domain size " << doSize << " elements, " << boxSize[0] << " by " << boxSize[1] << " by " << boxSize[2] << std::endl;

    auto domainObsPtr = domainObs.data().elements<int32_t>();
    aggregateDomain.assign(domainObsPtr, domainObsPtr + doSize);
    logger () << "Aggregate domain assigned" << std::endl;

    computeAndSendIncrementalBoundaries(aggregateDomain, previousSolid);
    aggregateDomain.clear();
    logger() << "Distributor done" << std::endl;

}


void Distributor::computeAndSendIncrementalBoundaries(const std::vector<int32_t>& domain, std::vector<int32_t>& previousSolid) {
        size_t len = domain.size();

        std::vector<int32_t> isSolid(len, 0); // 0 - fluid, 1 - solid, 2 - surface

        size_t n = 0;
        for (size_t i = 0; i < boxSize[X]; i++) {
            for (size_t j = 0; j < boxSize[Y]; j++) {
                for (size_t k = 0; k < boxSize[Z]; k++) {
                    if(domain[i*boxSize[Y]*boxSize[Z] + j*boxSize[Z] + k] != FLUID) {
//                        if(isSurfaceNode(i, j, k, domain)) {
//                            isSolid[n] = 2;
//                        }
//                        else {
//                            isSolid[n] = 1;
//                        }
                        isSolid[n] = 2;
                    }
                    n++;
                }
            }
        }
        logger() << "isSolid assigned" << std::endl;


        // int sizeIsSolid = (sizeof(isSolid) / sizeof(isSolid[0]));
        bool equal = true;

        if (isSolid.size() == previousSolid.size()) {
            for(size_t l = 0; l < isSolid.size(); ++l) {
                if (isSolid[l] != previousSolid[l]) {
                    equal = false;
                    break;
                }
            }
        } else {
            equal = false;
        }

        if (equal) {
            // Send empty message to signal: no change or further work.
            isSolid.clear();
            logger() << "Sending empty message: nothing to do" << std::endl;

        } else {
            // Store the changed solid mask
            previousSolid = isSolid;
            logger() << "Storing changed geometry" << std::endl;
        }
        auto solid_result = Data::grid(isSolid.data(), {boxSize[X], boxSize[Y], boxSize[Z]}, {"x", "y", "z"});

        instance.send("geometryOut1", Message(t_cur, solid_result));
        instance.send("geometryOut2", Message(t_cur, solid_result));
        logger() << "Distributor: data sent" << std::endl;
}

void runDistributor(int argc, char * argv[]) {
    Instance instance(argc, argv, {
            {Operator::F_INIT, {"domainIn"}},
            {Operator::O_F, {"geometryOut1",
                             "geometryOut2"}}
    });

    Distributor distributor(instance);

    while (instance.reuse_instance()) {
        distributor.execute();
    }
}

int main(int argc, char * argv[]) {
    runDistributor(argc, argv);
    return EXIT_SUCCESS;
}

