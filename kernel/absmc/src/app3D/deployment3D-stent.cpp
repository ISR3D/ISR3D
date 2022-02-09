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

    /// characteristic length scale
    slv.readCharLengthScales();

    /// load agents from file; includes vessel wall agents. Obstacles are loaded separately.
    slv.readDat3DFile("stage2." + slv.getBaseName() + ".dat");
    slv.initVerletNeighbourDetector();
    slv.readBondsFromFile("stage2." + slv.getBaseName() + "_nb.dat");

    // set agents at longitudinal (x) domain boundaries immobile.
    slv.setBoundaryMobility("deploy_boundary_mobility_x", "deploy_boundary_mobility_y", "deploy_boundary_mobility_z");

    slv.setDeploymentRules();

    slv.loadGeneratedStent();
    slv.loadGeneratedBalloon();

    slv.readIterationParameters("deploy_max_iter","deploy_vtk_iter","deploy_dat_iter");

    UnaryRadialDisplacementForce displForceStent = slv.createRadialStentDisplacementForce();

    TypeSpecificUnaryForce typeSpecificDisplacementForce(&slv.zeroForce);
    typeSpecificDisplacementForce.add(&displForceStent, tObstacle3D);

    slv.setForces(&slv.bForce, &slv.linAttrForce, &typeSpecificDisplacementForce);
    // smooth deployment: stent displacement is performed at each integration step

    slv.calculateCharForce();

    slv.initializeIntegrator();
    slv.setFixedIntervalController("deploy_iter_time");

    // vtp output parameters
    slv.readVtpParameters("deploy_vtp_scalars", "deploy_vtp_vectors");

    // main expansion cycle (until maxIter)
    slv.iterate();

    logger() << "Doing post-deployment equilibration..." << std::endl;

    /// Final iteration with a smaller eps to bring the system closer to the equilibrium
    /// reset forces, remove stent displacement force
    slv.setForces(&slv.bForce, &slv.linAttrForce, &slv.zeroForce);

    slv.readIterationParameters("deploy_final_eq_iter");
    slv.initializeIntegrator();
    slv.setFixedIntervalController("deploy_final_eq_time");

//    RungeKuttaIntegratorMPI<agent_t>
//        finalIntegrator(*uForce, *binForce, *bondForce, L, solverMaxDispl / 2, solverMaxDt * 10, nbDetector, integratorLogFile);

    slv.iterate();
    slv.removeBalloon();
    logger() << "Balloon removed, equilibrating..." << std::endl;
    slv.iterate();
    /// write the final configuration to file
    logger() << std::endl << "Deployment complete, writing output file" << std::endl;
    slv.writeAllOutputs("stage3.");

    return EXIT_SUCCESS;
}

