#include "absmc3D.h"

using namespace absmc;
const std::string solverName = "deployment";

/// main
int main(int argc, char *argv[])
{
    if (argc!=2) {
        std::cout << "Usage: " << solverName << " configDir/configName" << std::endl;
        return EXIT_FAILURE;
    }

    /// read config filename; output to the config directory
    const std::string configFileName(argv[1]);
    AgentSolver<3> slv(solverName);
    slv.loadConfig(configFileName);

    /// agent input files
    slv.readIOFileNames("deploy_input_file");
    slv.readCenterlineFromFile();
    slv.readCharLengthScales();

    /// load agents from file; includes vessel wall agents. Obstacles are loaded separately.
    slv.readDat3DFile("stage2." + slv.getBaseName() + ".dat");
    slv.initVerletNeighbourDetector();
    slv.readBondsFromFile("stage2." + slv.getBaseName() + "_nb.dat");

    slv.readWallDisplacements();
    slv.readStentWithDisplacements("deploy_stent_file", "deploy_stent_file_final", "displacements_stent");
    slv.setIterationParametersFromDisplacements("deploy_vtk_iter","deploy_dat_iter");

    /// set agents at longitudinal (x) domain boundaries immobile;
    slv.setBoundaryMobility("deploy_boundary_mobility_x", "deploy_boundary_mobility_y", "deploy_boundary_mobility_z");

    slv.detectSurface();
    slv.surfaceToMembranes();

    slv.setDeploymentRules();

    UnaryTrajectoryForce trajForce = slv.createWallTrajectoryForce();
    UnaryTrajectoryForce trajForceStent = slv.createStentTrajectoryForce();

    TypeSpecificUnaryForce typeSpecificDisplacementForce(&trajForce);
    typeSpecificDisplacementForce.add(&trajForceStent, tObstacle3D);

    slv.setForces(&slv.bForce, &slv.linAttrForce, &typeSpecificDisplacementForce);

    slv.calculateCharForce();


    slv.initializeIntegrator();
    slv.setFixedIntervalController("deploy_iter_time");

    // vtp output parameters
    slv.readVtpParameters("deploy_vtp_scalars", "deploy_vtp_vectors");

    slv.iterateWithTrajectories(trajForce, trajForceStent);

    /// write the final configuration to file
    std::cout << std::endl << "Deployment complete, writing output file" << std::endl;
    slv.writeAllOutputs("stage3.");

    return EXIT_SUCCESS;
}
