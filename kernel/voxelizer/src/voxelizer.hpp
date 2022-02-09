#ifndef VOXEL_HPP
#define VOXEL_HPP

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <limits>
#include <list>
#include <cstdio>
#include <array>
#include <stdexcept>

#include <libmuscle/libmuscle.hpp>
#include <ymmsl/ymmsl.hpp>

#include <util/logger.h>
#include <util/typeNames.hpp>
#include <mapping/sphereVoxelizer.hpp>

using libmuscle::Data;
using libmuscle::DataConstRef;
using libmuscle::Instance;
using libmuscle::Message;
using ymmsl::Operator;

using util::logger;

class VoxelizingMapper {
private:
    double length;
    double deltaX;
    double smooth_radius;
    int surface_iter;
    Instance& instance;
    const std::string outPath;

public:
    VoxelizingMapper(Instance& instance_):
    instance(instance_),
    outPath (instance.get_setting_as<std::string> ("muscle_data_out_dir")) {
        const std::string logName = instance.get_setting_as<std::string> ("log_file");
        if(!logName.empty())
            logger().setFile(outPath + logName);
    }

    void execute();

};

void VoxelizingMapper::execute() {
    static int iteration = 0;
    static std::vector<double> domainBoundingBox (6);
    static int detectorY = 0;
    static int detectorZ = 0;

    // A vector for all Sphere objects we receive this iteration
    std::vector<Sphere> spheres;

    // Read parameters
    deltaX = instance.get_setting_as<double>("flowdx");
    smooth_radius = instance.get_setting_as<double>("smooth_Radius");
    surface_iter = instance.get_setting_as<int64_t>("surface_iter");

    logger() << "Receiving cell position data from AB SMC submodel, iteration " << iteration << std::endl;
    // Receive message...
    auto msg = instance.receive("cellPositionsIn");
    logger() << "Received ..." << std::endl;
    auto spheresData = msg.data();
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

//    double x = spheresData[i][0][0].as<double>();
//    double y = spheresData[i][0][1].as<double>();
//    double z = spheresData[i][0][2].as<double>();
//    double r = spheresData[i][1].as<double>();
//    int typeID = spheresData[i][2].as<int64_t>(); // type of tissue; so far, it's smc or stent

    // A new Sphere object is made for each received item
    for (size_t i = 0; i < nAgents; ++i) {
        spheres.push_back(Sphere(posX[i],posY[i],posZ[i],rad[i],agentType[i]));
    }

    logger() << spheres.size() << " cells received." << std::endl;

    // make a voxel geometry
    logger() << "Voxelizing" << std::endl;
    // assuming the outer wall is static, we only have to initialize SurfaceDetector once at the 1st step
    if (iteration == 0) {
        double minX, minY, minZ;
        minX = minY = minZ = std::numeric_limits<double>::max();

        double maxX, maxY, maxZ;
        maxX = maxY = maxZ = std::numeric_limits<double>::lowest();

        //A bounding box is constructed which includes all cells.
        //It is assumed that the inlet and outlet are already perpendicular to X axis and there is no leakage
        // Look for max/min cell coordinates
        for (Sphere& s : spheres) {
            if (s.X <= minX) {
                minX = s.X;
            }
            if (s.X >= maxX) {
                maxX = s.X;
            }
            if (s.Y <= minY) {
                minY = s.Y;
            }
            if (s.Y >= maxY) {
                maxY = s.Y;
            }
            if (s.Z <= minZ) {
                minZ = s.Z;
            }
            if (s.Z >= maxZ) {
                maxZ = s.Z;
            }
        }

        //add padding around the vessel and adjust X-axis bounds
        maxX -= deltaX;
        minX += deltaX;
        maxY += 2 * deltaX;
        minY -= 2 * deltaX;
        maxZ += 2 * deltaX;
        minZ -= 2 * deltaX;
        logger() << "Domain: from (" << minX << "," << minY << "," << minZ << ") to (" << maxX << "," << maxY << "," << maxZ << ")." << std::endl;
        domainBoundingBox = {minX, minY, minZ, maxX, maxY, maxZ};

        // Look for the min and max of the lowest-X slice
        double lowestXminY = std::numeric_limits<double>::max();
        double lowestXminZ = std::numeric_limits<double>::max();
        double lowestXmaxY = std::numeric_limits<double>::lowest();
        double lowestXmaxZ = std::numeric_limits<double>::lowest();

        for (Sphere& s : spheres) {
            if (s.X <= minX + 2 * deltaX) {
                if (s.Y <= lowestXminY) {
                    lowestXminY = s.Y;
                }
                if (s.Y >= lowestXmaxY) {
                    lowestXmaxY = s.Y;
                }
                if (s.Z <= lowestXminZ) {
                    lowestXminZ = s.Z;
                }
                if (s.Z >= lowestXmaxZ) {
                    lowestXmaxZ = s.Z;
                }
            }
        }
        detectorY = static_cast<int>((lowestXmaxY + lowestXminY - 2 * minY) / (2. * deltaX));
        detectorZ = static_cast<int>((lowestXmaxZ + lowestXminZ - 2 * minZ) / (2. * deltaX));
    }

    // Initialize percolation detector with a node that's definitely inside the lumen
    PercolationDetector surfaceDetector (int_point(0, detectorY, detectorZ));

    logger() << "Creating new sphereVox" << std::endl;
    SphereVoxelizer sphereVox(deltaX, domainBoundingBox, smooth_radius, surfaceDetector);

    logger() << "Set values of sphereVox" << std::endl;

    sphereVox.setShapes(spheres);
    sphereVox.voxelize();
    sphereVox.percolateOutside();

    //set stent lattice points from "solid" to "strut"
    std::vector<int32_t> aggregateDomain;
    aggregateDomain = sphereVox.getSmoothDomain();
    logger() << "Smooth domain calculated" << std::endl;

    double t_cur = msg.timestamp();

    auto box_size_result = Data::grid(domainBoundingBox.data(), {domainBoundingBox.size()}, {"xyz"});
    instance.send("domainBoxOut", Message(t_cur, box_size_result));

    //Count nodes by type, init all types with zeroes
    std::array <size_t, 256> numNodes;
    numNodes.fill(0);
    for (auto val : aggregateDomain) {
        if (val >= 256) {
            throw std::runtime_error("More than 255 tissue types are not supported");
        }
        numNodes[val]++;
    }

    logger() << "Total nodes: " << (aggregateDomain.size()) <<
            ", fluid nodes: " << numNodes[FLUID] << ", solid: " << numNodes[SOLID] << ", strut: " << numNodes[STRUT] << ", static nodes: " << numNodes[STATIC] << std::endl;

    /// no char type in Muscle3 (yet)
    //std::vector<int32_t> aggregateDomainForSend (aggregateDomain.begin(), aggregateDomain.end());

    auto aD_result = Data::grid(aggregateDomain.data(), {sphereVox.sizeX, sphereVox.sizeY, sphereVox.sizeZ}, {"x", "y", "z"});
    //auto aD_result = Data::grid(testVector.data(), {100, 100, 100}, {"x", "y", "z"});

    instance.send("domainOut1", Message(t_cur, aD_result));
    instance.send("domainOut2", Message(t_cur, aD_result));
    logger() << "DomainEntrances sent" << std::endl;

    std::vector<double> lumenArea = sphereVox.getLumenArea();

    char str[7];
    snprintf (str, 7, "%06d", iteration);

    if((surface_iter != 0) && (iteration % surface_iter == 0)) {
        sphereVox.writeSurfaceToFile(outPath + "/surface_points_" + std::string(str) + ".xyz");
    }

    std::string lumenFilename = outPath + "/lumen_area_" + std::string(str) + ".csv";
    try{
       std::ofstream file;
       file.open(lumenFilename);
       std::string dlmtr = "";
       for (auto area : lumenArea) {
           file << dlmtr << area;
           dlmtr = '\t';
       }
       file << std::endl;
       file.close();
    } catch (std::runtime_error& e) {
        logger() << "RUNTIME ERROR, CANNOT WRITE LUMEN FILE: " << e.what() << std::endl;
    }

    instance.send("cellsOut", Message(t_cur, spheresData));

    std::vector<int_point> surf = sphereVox.getSurface();
    auto surf_result = Data::nils(surf.size());
    for (std::size_t i = 0; i < surf.size(); ++i) {
        surf_result[i] = Data::list(surf[i][0], surf[i][1], surf[i][2]);
    }
    instance.send("surfacePositionOut", Message(t_cur, surf_result));
    logger() << "Sent surface positions to the collector" << std::endl;

    //sphereVox->dispose();
    iteration++;
    logger() << "Out of cell mapper loop" << std::endl;
}

#endif
