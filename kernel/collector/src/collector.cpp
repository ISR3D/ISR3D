#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <vector>
#include <stdexcept>

#include <libmuscle/libmuscle.hpp>
#include <ymmsl/ymmsl.hpp>

#include <util/logger.h>
#include <util/typeNames.hpp>

#include <mapping/sphere.hpp>
#include <mapping/sphereVoxelizer.hpp>
#include <mapping/intPoint.hpp>
#include <mapping/surfaceDetector.hpp>

using libmuscle::Data;
using libmuscle::DataConstRef;
using libmuscle::Instance;
using libmuscle::Message;
using ymmsl::Operator;

using util::logger;

class Collector {
private:
    double deltaX;
    double t_cur;
    Instance& instance;
    const std::string outPath;
    SphereVoxelizer* spherevox = nullptr;

    std::vector<double> previousShearPalabosData;
    std::vector<int32_t> previousSolid;
    std::vector<double> drug;
    std::vector<double> neointima;
    std::vector<double> cellShearpalabos;
    std::vector<double> domainBoundingBox; //{minX, minY, minZ, maxX, maxY, maxZ}
    std::vector<size_t> boxSize; // number of voxels along X, Y, Z

public:
    void reset();
    void execute();

    Collector(Instance& instance_):
    instance(instance_),
    outPath (instance.get_setting_as<std::string> ("muscle_data_out_dir")) {
        const std::string logName = instance.get_setting_as<std::string> ("log_file");
        if(!logName.empty())
            logger().setFile(outPath + logName);

        deltaX = instance.get_setting_as<double>("flowdx");
        logger() << "Collector initialised" << std::endl;
    }

};

void Collector::execute() {
    auto d_msg = instance.receive("domainBoxIn");
    if (spherevox == nullptr) {
        auto domainPtr = d_msg.data().elements<double>();
        domainBoundingBox.assign(domainPtr, domainPtr + d_msg.data().size());

        logger() << "Creating voxelizer" << std::endl;
        NoSurface surfaceDetector;
        spherevox = new SphereVoxelizer(deltaX, domainBoundingBox, 0., surfaceDetector);
    }

    // to receive cell lists from SMC
    logger() << "Receiving cell position" << std::endl;
    auto c_msg = instance.receive("cellsIn");

    logger() << "Received ..." << std::endl;
    auto spheresData = c_msg.data();
    // ...and parse all received data fields

    size_t nAgents = spheresData[0].size();
    std::vector<double> posX(nAgents), posY(nAgents), posZ(nAgents), rad(nAgents);
    std::vector<int64_t> agentType(nAgents);

    auto posXPtr = spheresData[0].elements<double>();
    posX.assign(posXPtr, posXPtr + nAgents);

    auto posYPtr = spheresData[1].elements<double>();
    posY.assign(posYPtr, posYPtr + nAgents);

    auto posZPtr = spheresData[2].elements<double>();
    posZ.assign(posZPtr, posZPtr + nAgents);

    auto radPtr = spheresData[3].elements<double>();
    rad.assign(radPtr, radPtr + nAgents);

    auto typePtr = spheresData[4].elements<int64_t>();
    agentType.assign(typePtr, typePtr + nAgents);

    // A new Sphere object is made for each received item
    std::vector<Sphere> spheres;
    for (size_t i = 0; i < nAgents; ++i) {
        spheres.push_back(Sphere(posX[i],posY[i],posZ[i],rad[i],agentType[i]));
    }

    logger() << spheres.size() << " cells received." << std::endl;
    spherevox->setShapes(spheres);

    //receive mesh domain
    //similar to distributor, can be moved to a function
    auto domainObs = instance.receive("domainIn");
    t_cur = domainObs.timestamp();
    // nextTime = domainObs.getNextTimestamp();
    logger() << "Domain received for time " << t_cur << std::endl;

    size_t doSize = domainObs.data().size();
    boxSize = domainObs.data().shape();
    logger() << "Domain size " << doSize << " elements, " << boxSize[0] << " by " << boxSize[1] << " by " << boxSize[2] << std::endl;

    std::vector<int32_t> aggregateDomain;
    auto domainObsPtr = domainObs.data().elements<int32_t>();
    aggregateDomain.assign(domainObsPtr, domainObsPtr + doSize);
    /// convert to int32_t for voxelizer
    //std::vector<int32_t> aggregateDomainForVox (aggregateDomain.begin(), aggregateDomain.end());
    spherevox->setMaskedDomain(aggregateDomain);
    logger () << "Sphere voxelizer domain assigned" << std::endl;

    //can't receive references
    //auto cp_msg = instance.receive("cellPositionExit");
    //spherevox->setShapeReference(cp.data());
    //rebuild voxelization locally instead
    spherevox->voxelizeNoSurface();

    auto sp_msg = instance.receive("surfacePositionsIn");
    auto spData = sp_msg.data();
    std::vector<int_point> surfacePoints;
    for (size_t i = 0; i < spData.size(); ++i) {
        surfacePoints.push_back(int_point(spData[i][0].as<int32_t>(),
                                          spData[i][1].as<int32_t>(),
                                          spData[i][2].as<int32_t>()));
    }

    logger() << "Receiving BF data" << std::endl;
    spherevox->setSurface(surfacePoints);

    //log("Mapping neointima distance to strut");
    //neointima = spherevox->mapDomain(neointimaDistances(domain, size), 0.1);

    //receive drug from DD
    //log("Receiving drug concentration.");
    //Observation<double[][][]> concObs = concentrationExit.receiveObservation();
    //double[][][] drugConcentration = concObs.getData();

    //log("Mapping drug concentration");
    //drug = spherevox->mapDomain(drugConcentration, 0.0);

    //receive surface values from BF
    logger() << "Receiving BF data" << std::endl;

    auto coord_msg = instance.receive("geometryIn");
    std::vector<int32_t> isSolid;
    auto solPtr = coord_msg.data().elements<int32_t>();
    isSolid.assign(solPtr, solPtr + coord_msg.data().size());

    auto palabos_msg = instance.receive("shearStressIn");
    std::vector<double> shearpalabosData;
    auto shearPtr = palabos_msg.data().elements<double>();
    shearpalabosData.assign(shearPtr, shearPtr + palabos_msg.data().size());

    if (isSolid.size() != 0) {
        logger() << "New BF results" << std::endl;
        // Store the changed values
        previousSolid = isSolid;
        previousShearPalabosData = shearpalabosData;
    } else {
        // If no change happened, use previous values
        logger() << "No change in BF results, using previous results" << std::endl;
        isSolid = previousSolid;
        shearpalabosData = previousShearPalabosData;
    }

    if (isSolid.size() != shearpalabosData.size()) {
        throw std::runtime_error(std::string("Size of data mask (") + std::to_string(isSolid.size()) + ") is not equal to the size of the data (" + std::to_string(shearpalabosData.size()) + ")");
    }
    if (isSolid.size() != boxSize[X]*boxSize[Y]*boxSize[Z]) {
        throw std::runtime_error(std::string("Size of data mask (") + std::to_string(isSolid.size()) + ") is not equal to the size of the grid (" + std::to_string((boxSize[X]*boxSize[Y]*boxSize[Z])) + ")");
    }

    logger() << "Mapping wall shear stress" << std::endl;
    /// convert to int32_t for voxelizer
    std::vector<int32_t> isSolidForVox (isSolid.begin(), isSolid.end());
    cellShearpalabos = spherevox->mapSurfaceAverage(shearpalabosData, isSolidForVox);

//        log("Drug size : " + drug.length);
//        log("Neointima size : " + neointima.length);
    logger() << "WSS  size : " << cellShearpalabos.size() << std::endl;

    //send mapped data to SMC
    logger() << "Sending data to SMC." << std::endl;
    auto cell_result = Data::grid(cellShearpalabos.data(), {cellShearpalabos.size()});
    instance.send("cellMaxStressOut", Message(t_cur, cell_result));
//        NeoEntrance.send(neointima, time, nextTime);
    logger() << "Data to SMC sent." << std::endl;

    reset();

}

void Collector::reset() {
    t_cur = 0.;
    // nextTime = NULL;
    cellShearpalabos.clear();
//        neointima = null;
    spherevox->dispose();
}

void runCollector(int argc, char * argv[]) {
    Instance instance(argc, argv, {
            {Operator::F_INIT, {"cellsIn",
                                "domainBoxIn",
                                "domainIn",
                                "surfacePositionsIn",
                                "shearStressIn",
                                "geometryIn"}},
            {Operator::O_F, {"cellMaxStressOut"}}});

    Collector collector(instance);

    while (instance.reuse_instance()) {
        collector.execute();
    }
}

int main(int argc, char * argv[]) {
    runCollector(argc, argv);
    return EXIT_SUCCESS;
}
