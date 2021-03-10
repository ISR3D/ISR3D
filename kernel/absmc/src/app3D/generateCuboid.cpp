#include <iostream>
#include <vector>
#include <string>

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

void cuboidPointsToCells(std::vector<point_t> const& points,
                                 AgentContainer<3> & agentContainer,
                                 AgentFactory<3> const& agentFactory,
                                 double r,
                                 double dSplit, double dSplitout)
{
    const double dSplitSqr = dSplit*dSplit;
    const double dSplitoutSqr = dSplitout*dSplitout;

    const size_t nPoints = points.size();
    for (size_t iPoint=0; iPoint<nPoints; iPoint++) {

        point_t const& p = points[iPoint];
        const double dSqr = p[1]*p[1]  + p[2]*p[2];

        const AgentTypeId typeId =   (p[2] > 0.5-r-0.00001) ? tEEL3D
        : tSMC3D ;


        //        if (p[0] > r+0.00001) {
        AgentBase<3>* newAgent = agentFactory.create(typeId, p, p, r);
        //        }


        if (newAgent) {
            //            newAgent ->setGhostIndex(agentContainer.size());
            agentContainer.add(newAgent);
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc!=2) {
        cout << "Usage: generateCuboid configDir" << endl;
        return EXIT_FAILURE;
    }

    const std::string configDir(argv[1]);
    const std::string configFileName  = createFileName(configDir, "absmc", "cfg");
    config().readFile(configFileName);

    // domain dimensions
    const double lX = config().getValue<double>("lX");
    const double lY = config().getValue<double>("lY");
    const double lZ = config().getValue<double>("lZ");

    // morse potential parameters (dimensionless)
    const double U1 = config().getValue<double>("U1"); // 5.0;
    const double r1 = config().getValue<double>("r1"); // 1.0;
    const double U2 = config().getValue<double>("U2"); // 5.0;
    const double r2 = config().getValue<double>("r2"); // 5.0;
    // equilibrium distance (dimensionless)
    const double r0 = r1*r2/(r2-r1)*log(U1/U2*r2/r1);
    // characteristic length scale
    const double L = config().getValue<double>("smc_mean_rad");

    const double meanShift  = config().getValue<double>("generate_sigma_shift");
    const double sigmaShift = config().getValue<double>("generate_sigma_shift");

    std::vector<point_t> pointsDRH;
    GeometryGenerator3D::generateRegularCuboidDRH(pointsDRH, r0*L, lX, lY, lZ);
    GeometryGenerator3D::shufflePositions(pointsDRH, meanShift, sigmaShift);

    AgentContainer<3> agentsDRH(pointsDRH.size());
    pointsToCells<3>(pointsDRH, agentsDRH, AgentFactory3D(), tSMC3D, 0.5*r0*L);
    agentsDRH.writeFile(createFileName(configDir, "cuboid_drh", "dat") );
    VtkWriter<agent_t>::writeVtkPolyData(agentsDRH.getAgentVector(), createFileName(configDir, "cuboid_drh", "vtp"), "typeId", "");


    std::vector<point_t> pointsJez;
    std::vector<index_vec_t> nb;
    GeometryGenerator3D::generateRegularCuboidJez(pointsJez, r0*L, lX, lY, lZ,&nb);
    GeometryGenerator3D::shufflePositions(pointsJez, meanShift, sigmaShift);

    AgentContainer<3> agentsJez(pointsJez.size());
 //   pointsToCells<3>(pointsJez, agentsJez, AgentFactory3D(), tSMC3D, 0.5*r0*L);
    cuboidPointsToCells(pointsJez, agentsJez, AgentFactory3D(), 0.5*r0*L, 0.0, 0.0);
    agentsJez.writeFile(createFileName(configDir, "cuboid_jez", "dat") );
    VtkWriter<agent_t>::writeVtkPolyData(agentsJez.getAgentVector(), createFileName(configDir, "cuboid_jez", "vtp") );

    NBWriter<3>::writeFile(nb, createFileName(configDir, "smooth_cylinder_nb", "dat") );
    VtkWriter<agent_t>::writeVtkPolyData(agentsJez.getAgentVector(), createFileName(configDir, "smooth_cylinder", "vtp"), "typeId", "");

    // create simple stent consisting of num_struts straight wires
    cout << "Generating \" strut..." << endl;
    //const int numStruts = config().getValue<int>("num_struts");

    std::vector<point_t> strutPoints;
    GeometryGenerator3D::generateSimpleStrut(strutPoints, r0*L, lX, -r0*L);

    AgentContainer<3> strutAgents(strutPoints.size() );
    pointsToCells<3>(strutPoints, strutAgents, AgentFactory3D(), tObstacle3D, 0.5*r0*L);

    strutAgents.writeFile(createFileName(configDir, "simple_stent", "dat"), tObstacle3D);
    VtkWriter<agent_t>::writeVtkPolyData(strutAgents.getAgentVector(), createFileName(configDir, "simple_stent", "vtp"), "typeId", "");

    return EXIT_SUCCESS;
}
