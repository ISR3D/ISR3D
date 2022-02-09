#include <iostream>

#include "absmc3D.h"

using namespace absmc;

const std::string solverName = "equilibration";

/// main
int main(int argc, char *argv[])
{
    if (argc!=2) {
        std::cout << "Usage: " << solverName << " configDir/configName" << std::endl;
        return EXIT_FAILURE;
    }

    const std::string configFileName(argv[1]);

    AgentSolver<3> slv(solverName);
    slv.loadConfig(configFileName);

    // agent input file
    slv.readIOFileNames("eq_input_file");
    slv.readCharLengthScales();

    // load agents from file includes smc+iel+obstacle agents, specify prefix here
    slv.readDat3DFile("stage1." + slv.getBaseName() + ".dat");

    /// initialise Verlet neighbourhoods
    slv.initVerletNeighbourDetector();

    /// Read bonds
    slv.readBondsFromFile("stage1." + slv.getBaseName() + "_nb.dat");

    // set agents at longitudinal (x) domain boundaries immobile.
    slv.setBoundaryMobility("eq_boundary_mobility_x", "eq_boundary_mobility_y", "eq_boundary_mobility_z");

    // randomly displace agents; only the specified fraction will be displaced,
    // the remainder is left in position and set immobile
    slv.shuffleAgentsByConfig();

    slv.setForces(&slv.bForce, &slv.linAttrForce, &slv.zeroForce);
    slv.calculateCharForce();
    slv.initializeIntegrator();
    slv.readIterationParameters("eq_max_iter","eq_vtk_iter","eq_dat_iter");
    slv.setForceResidualMaxNormController("eq_convergence_level");

    // vtp output parameters
    slv.readVtpParameters("eq_vtp_scalars", "eq_vtp_vectors");

    logger() << "Starting main loop..." << std::endl;
    // main eqilibration cycle
    slv.iterate();

    // create bonds from the detected neighbourhoods
    slv.makeBondsFromNeighbours();
    // set initial positions of cells to current positions
    slv.resetPos0();
    // write the final configuration to file
    slv.writeAllOutputs("stage2.");

    return EXIT_SUCCESS;
}
