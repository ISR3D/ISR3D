#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

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


/// Create agents from the given vector of points using the given factory, assuming constant radius.
/// For points with distance from the cylindrical axis smaller than dSplit, an IEL3D agent is created, otherwise a SMC3D one.
void smoothCylinderPointsToCells(std::vector<point_t> const& points,
                                 AgentContainer<3> & agentContainer,
                                 AgentFactory<3> const& agentFactory,
                                 double r, double sigma,
                                 double dSplit, double dSplitout)
{
    const double dSplitSqr = dSplit*dSplit;
    const double dSplitoutSqr = dSplitout*dSplitout;

    util::Random random;

    const size_t nPoints = points.size();
    for (size_t iPoint=0; iPoint<nPoints; iPoint++) {

        point_t const& p = points[iPoint];
        const double dSqr = p[1]*p[1]  + p[2]*p[2];
        double curR = random.getNextGauss(r, sigma);
        const AgentTypeId typeId =  (dSqr < dSplitSqr) ? tIEL3D
                                    : (dSqr > dSplitoutSqr) ? tEEL3D
                                    : tSMC3D ;

        AgentBase<3>* newAgent = agentFactory.create(typeId, p, p, curR);

        if (newAgent) {
            agentContainer.add(newAgent);
        }
    }
}

/// bend the vessel into two semicircles with a straight segment in the middle
void bend(AgentContainer<3> & agents, const double curvatureRadius, const double lX, const double straightLength = 0.0) {
    point_t centerLow ((lX - straightLength)/2, -curvatureRadius, 0); //center for low-X semicircle
    point_t centerHigh ((lX + straightLength)/2, -curvatureRadius, 0); //center for high-X semicircle
    size_t nAgents = agents.count();
    for (size_t iAgent = 0; iAgent < nAgents; iAgent++) {
        point_t agentPos = agents.getAgentVector()[iAgent]->getPos();
        double agentR = agents.getAgentVector()[iAgent]->getR();
        if ((agentPos[0] < centerLow[0])||(agentPos[0] > centerHigh[0])) {
            agents.getAgentVector()[iAgent]->setR(agentR * (curvatureRadius + agentPos[1]) / curvatureRadius);
            point_t center;
            if (agentPos[0] < centerLow[0]) {
                center = centerLow;
            }
            if (agentPos[0] > centerHigh[0]) {
                center = centerHigh;
            }
            const double phi = (agentPos[0] - center[0]) / curvatureRadius;
            const double xPrime = (agentPos[1]+curvatureRadius) * sin(phi) + center[0];
            const double yPrime = (agentPos[1]+curvatureRadius) * cos(phi) + center[1];
            agentPos[0] = xPrime;
            agentPos[1] = yPrime;
            agents.getAgentVector()[iAgent]->setPos(agentPos);
        }
    }
}

/// Create "smooth" version of cylindrical shell geometry by wrapping a cuboid around a cylinder.
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
    const double innerRadiusConfig = config().getValue<double>("inner_radius");
    const double outerRadius = config().getValue<double>("outer_radius");

    // equilibrium distance (dimensionless)
    const double r0 = 2;
    // characteristic length scale
    const double L = config().getValue<double>("smc_mean_rad");
    // Radius sigma (standard deviation)
    const double sigma = config().getValue<double>("smc_rad_sigma");
    const std::string outputName = config().getValue<std::string>("generate_out_name");
    const std::string smoothCylinderName = parseFileName(outputName).baseName;

    std::vector<point_t> points;
    std::vector<point_t> balloonPoints;
    std::vector<index_vec_t> nb;

    cout << "Generating \"smooth\" cylinder..." << endl;
    /// ensure that wall thickness is a multiple of cell spacing
    const double dZ = r0*L / 2;
    double nLayers = floor((outerRadius - innerRadiusConfig) / dZ);
    if ((int)nLayers % 2 == 1)
        nLayers += 1.0;
    const double innerRadius = outerRadius - nLayers * dZ;

    GeometryGenerator3D::generateSmoothCylinder(points, r0*L, lX, innerRadius, outerRadius, &nb);

    std::cout << std::endl;

    const double meanShift  = config().getValue<double>("generate_sigma_shift");
    const double sigmaShift = config().getValue<double>("generate_sigma_shift");
    GeometryGenerator3D::shufflePositions(points, meanShift, sigmaShift);

    AgentContainer<3> agents(points.size() );
    /// 0.1 gaps to avoid comparing almost equal doubles
    smoothCylinderPointsToCells(points, agents, AgentFactory3D(), L, sigma, innerRadius+0.6*r0*L, outerRadius-0.6*r0*L); //outerRadius-0.5*r0*L);

    /// bend the vessel
    const double curvatureRadius = config().getValue<double>("curvature_radius");
    const double straightLength = config().getValue<double>("straight_segment");
    bend(agents, curvatureRadius, lX, straightLength);

    cout << "Writing files for a simple artery..." << endl;
    VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(), createFileName(configDir, "stage1."+smoothCylinderName, "vtp"), "typeId", "");
    agents.writeFile(createFileName(configDir, "stage1."+smoothCylinderName, "dat") );
    NBWriter<3>::writeFile(nb, createFileName(configDir, "stage1."+smoothCylinderName + "_nb", "dat") );

    cout << agents.count() << " agents generated." << endl;

    // create simple stent consisting of num_struts straight wires

    const int numStruts = config().getValue<int>("num_struts");
    const int stentLength = config().getValue<double>("stent_length");
    const int helixPitch = config().getValue<double>("helix_pitch");
    const double stentCurvatureRadius = config().getValue<double>("stent_curvature_radius");

    std::vector<point_t> strutPoints;
    GeometryGenerator3D::generateSimpleStent(strutPoints, L, stentLength, innerRadius-r0*L, numStruts, helixPitch);

    cout << "Generating a simple stent; " << strutPoints.size() << " points in the geometry."<< endl;
    AgentContainer<3> strutAgents(strutPoints.size() );
    pointsToCells<3>(strutPoints, strutAgents, AgentFactory3D(), tObstacle3D, 0.5*r0*L);
    bend(strutAgents, stentCurvatureRadius, stentLength);

    strutAgents.writeFile(createFileName(configDir, "simple_stent", "dat"), tObstacle3D);
    VtkWriter<agent_t>::writeVtkPolyData(strutAgents.getAgentVector(), createFileName(configDir, "simple_stent", "vtp"), "typeId", "");

    /// generate a simple "balloon"
    GeometryGenerator3D::generateSmoothCylinder(balloonPoints, r0*L, stentLength, innerRadius - 2 * r0*L - 0.050, innerRadius - r0*L - 0.050);
    /// Alternative: use a specific balloon generator:
    /// SimpleBalloon3D balloon(balloonAgents, yOffset, 0.0, balloonXMin, balloonXMax, L * 2);

    cout << "Generating a simple balloon; " << balloonPoints.size() << " points in the geometry." << endl;
    AgentContainer<3> balloonAgents(balloonPoints.size() );
    pointsToCells<3>(balloonPoints, balloonAgents, AgentFactory3D(), tObstacle3D, r0*L);
    bend(balloonAgents, stentCurvatureRadius, stentLength);

    balloonAgents.writeFile(createFileName(configDir, "simple_balloon", "dat"), tObstacle3D);
    VtkWriter<agent_t>::writeVtkPolyData(balloonAgents.getAgentVector(), createFileName(configDir, "simple_balloon", "vtp"), "typeId", "");

    /// Create a centerline file
    /// Same curvature as stent and balloon, used for deployment
    cout << "Generating a centerline file." << endl;
    std::vector<point_t> centerlinePoints;
    double lineStep = 0.05;
    for (double xi = 0.; xi <= lX; xi += lineStep) {
        if(stentCurvatureRadius != 0) {
            point_t center (lX/2, -stentCurvatureRadius, 0);
            const double phi = (xi - center[0]) / stentCurvatureRadius;
            const double xPrime = (stentCurvatureRadius) * sin(phi) + center[0];
            const double yPrime = (stentCurvatureRadius) * cos(phi) + center[1];
            centerlinePoints.push_back(point_t(xPrime, yPrime, 0.0));
        }
        else {
            centerlinePoints.push_back(point_t(xi, 0.0, 0.0));
        }
    }
    Centerline<3> centerline;
    centerline.initializeFromVector(centerlinePoints);
    const std::string centerlineFileName = config().getValue<std::string>("centerline_file");
    centerline.writeToFile(createFileName(configDir,centerlineFileName));

    cout << "Generating geometries done." << endl;

    return EXIT_SUCCESS;
}
