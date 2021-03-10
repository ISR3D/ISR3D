#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "absmc3D.h"

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

/// Create "simple" version of cylindrical shell geometry by cutting a shell out of a cuboid.
int main(int argc, char *argv[])
{
    if (argc!=2) {
        cout << "Usage: generateSimpleCylinder configDir" << endl;
        return EXIT_FAILURE;
    }

    const std::string configDir(argv[1]);
    const std::string configFileName  = createFileName(configDir, "absmc", "cfg");
    config().readFile(configFileName);

    // domain size (real dimensions)
    const double lX = config().getValue<double>("lX");
    const double innerRadius = config().getValue<double>("inner_radius");
    const double outerRadius = config().getValue<double>("outer_radius");

    // morse potential parameters (dimensionless)
    const double U1 = config().getValue<double>("U1");
    const double r1 = config().getValue<double>("r1");
    const double U2 = config().getValue<double>("U2");
    const double r2 = config().getValue<double>("r2");
    // equilibrium distance (dimensionless)
    const double r0 = r1*r2/(r2-r1)*log(U1/U2*r2/r1);
    // characteristic length scale
    const double L = config().getValue<double>("smc_mean_rad");

    std::vector<point_t> points;
    std::vector<index_vec_t> nb;

    cout << "Generating \"simple\" cylinder..." << endl;
    GeometryGenerator3D::generateSimpleCylinder(points, r0*L, lX, innerRadius, outerRadius, &nb);

    const double meanShift  = config().getValue<double>("generate_sigma_shift");
    const double sigmaShift = config().getValue<double>("generate_sigma_shift");
    GeometryGenerator3D::shufflePositions(points, meanShift, sigmaShift);

    cout << "Writing files for \"simple\" cylinder..." << endl;
    AgentContainer<3> agents(points.size() );
    pointsToCells<3>(points, agents, AgentFactory3D(), tSMC3D, 0.5*r0*L);
    agents.writeFile(createFileName(configDir, "simple_cylinder", "dat") );
    NBWriter<3>::writeFile(nb, createFileName(configDir, "simple_cylinder_nb", "dat") );
    VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(), createFileName(configDir, "simple_cylinder", "vtp"), "typeId", "");

    // create simple stent consisting of num_struts straight wires
    cout << "Generating \"simple\" stent..." << endl;
    const int numStruts = config().getValue<int>("num_struts");

    std::vector<point_t> strutPoints;
    GeometryGenerator3D::generateSimpleStent(strutPoints, r0*L, lX, innerRadius-r0*L, numStruts);

    AgentContainer<3> strutAgents(strutPoints.size() );
    pointsToCells<3>(strutPoints, strutAgents, AgentFactory3D(), tObstacle3D, 0.5*r0*L);

    strutAgents.writeFile(createFileName(configDir, "simple_stent", "dat"), tObstacle3D);
    VtkWriter<agent_t>::writeVtkPolyData(strutAgents.getAgentVector(), createFileName(configDir, "simple_stent", "vtp"), "typeId", "");

    return EXIT_SUCCESS;
}
