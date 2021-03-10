#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "absmc3D.h"
#include <util/random.h>

using namespace absmc;
using namespace absmc::graphics;

using util::createFileName;
using util::parseFileName;
using util::config;

using std::cout;
using std::cerr;
using std::endl;

typedef Point<3, double> point_t;
typedef AgentBase<3> agent_t;
typedef EuclideanMetrics<3, double> metrics_t;
typedef GeometryGenerator3D::index_vec_t index_vec_t;


/// Create agents from the given vector of points using the given factory
/// The type of agents is set as a parameter
/// Agent size can be randomized by setting non-zero sigma to generate uneven surfaces.
void pointsToCells(std::vector<point_t> const& points,
                                 AgentContainer<3> & agentContainer,
                                 AgentFactory<3> const& agentFactory, AgentTypeId typeId,
                                 double r, double sigma)
{

    util::Random random;

    const size_t nPoints = points.size();
    for (size_t iPoint=0; iPoint<nPoints; iPoint++) {

        point_t const& p = points[iPoint];
        double curR = random.getNextGauss(r, sigma);

        AgentBase<3>* newAgent = agentFactory.create(typeId, p, p, curR);

        if (newAgent) {
            agentContainer.add(newAgent);
        }
    }

}

/// Create a rectangular tunnel with ridges perpendicular to flow on the bottom surface.
int main(int argc, char *argv[])
{
    if (argc!=2) {
        cout << "Usage: generateSmoothCylinder configDir/configName" << endl;
        return EXIT_FAILURE;
    }

    const std::string configFileName(argv[1]);
    const std::string configDir = parseFileName(configFileName).directory;
    config().readFile(configFileName);

    // domain size (real dimensions)
    const double lX = config().getValue<double>("lX");
    const double lY = config().getValue<double>("lY");
    const double lZ = config().getValue<double>("lZ");
    const double ridge_height = config().getValue<double>("ridge_height");
    const double ridge_width = config().getValue<double>("ridge_width");

    const double ridge_firstX = config().getValue<double>("ridge_firstX");
    const double ridge_spacing = config().getValue<double>("ridge_spacing");
    const int ridge_number = config().getValue<int>("ridge_number");

    // equilibrium distance (dimensionless)
    const double r0 = 2;
    // characteristic length scale
    const double obstacle_rad = config().getValue<double>("obstacle_rad");
    //Radius sigma (standard deviation)
    const double sigma = 0.0;
    const std::string outputName = config().getValue<std::string>("generate_out_name");
    const std::string baseName = parseFileName(outputName).baseName;

    std::vector<point_t> points;
    std::vector<index_vec_t> nb;

    double wall_thickness = 0.015;
    cout << "Generating ridged rectangular channel." << endl;
    GeometryGenerator3D::generateRidgedChannel (points, obstacle_rad, lX, lY, lZ, ridge_height, ridge_width,
                                                ridge_firstX, ridge_spacing, ridge_number, wall_thickness);

    const double meanShift  = config().getValue<double>("generate_sigma_shift");
    const double sigmaShift = config().getValue<double>("generate_sigma_shift");
    GeometryGenerator3D::shufflePositions(points, meanShift, sigmaShift);

    const double endo_lowX = config().getValue<double>("endo_lowX");
    const double endo_highX = config().getValue<double>("endo_highX");
    const double endo_lowZ = config().getValue<double>("endo_lowZ");
    const double endo_highZ = config().getValue<double>("endo_highZ");
    const double endo_rad = config().getValue<double>("endo_rad");
    const double endo_spacing = config().getValue<double>("endo_spacing");
    std::vector<point_t> endo_points;
    GeometryGenerator3D::generateHexagonalMonolayer(endo_points, endo_spacing, -lY/2 + endo_rad,
                            endo_lowX, endo_lowZ, endo_highX, endo_highZ);

    cout << "Writing files for ridged rectangular channel." << endl;
    AgentContainer<3> agents(points.size() + endo_points.size() );
    pointsToCells(points, agents, AgentFactory3D(), tObstacle3D, obstacle_rad, sigma);
    pointsToCells(endo_points, agents, AgentFactory3D(), tEndo3D, endo_rad, sigma);
    agents.writeFile(createFileName(configDir, baseName, "dat") );
    NBWriter<3>::writeFile(nb, createFileName(configDir, baseName + "_nb", "dat") );
    VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(), createFileName(configDir, baseName, "vtp"), "typeId", "");

    cout << agents.count() << " agents loaded." << endl;


    return EXIT_SUCCESS;
}
