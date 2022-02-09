#include "absmc3D.h"
#include <vector>

using namespace absmc;
using util::config;

const std::string solverName = "uniaxialStrainTest";

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
    slv.readIOFileNames("strain_input_file");
    /// characteristic length scale
    slv.readCharLengthScales();
    slv.readIterationParameters("strain_max_iter");

    /// load agents from file; includes vessel wall agents. Obstacles are loaded separately.
    slv.readDat3DFile("stage2." + slv.getBaseName() + ".dat");
    slv.initVerletNeighbourDetector();
    slv.readBondsFromFile("stage2." + slv.getBaseName() + "_nb.dat");

    // set agents at longitudinal (x) domain boundaries immobile.
    slv.setBoundaryMobility("strain_boundary_mobility_x", "strain_boundary_mobility_y", "strain_boundary_mobility_z", "strain_boundary_immobility_range");

    slv.setDeploymentRules();

    Poly6thPowerAttractionForce3D bondForce = slv.createPolynomialForce("strain_force_c1","strain_force_c2","strain_force_c3","strain_force_c4","strain_force_c5","strain_force_c6");
    slv.setForces(&slv.bForce, &bondForce, &slv.zeroForce);

    slv.calculateCharForce();

    slv.initializeIntegrator();
    slv.setForceResidualMaxNormController("strain_convergence_level");

    // vtp output parameters
    slv.readVtpParameters("strain_vtp_scalars", "strain_vtp_vectors");

    const double strain_max_value = config().getValue<double>("strain_max_value");
    const double strain_deform_step = config().getValue<double>("strain_deform_step");
    std::vector<double> stressValues;

    for(double cur_strain = 0.0; cur_strain <= strain_max_value; cur_strain += strain_deform_step) {
        std::cout << "Strain: " << cur_strain << std::endl;
        slv.applyEngineeringStrain(cur_strain);
        // main expansion cycle (until maxIter)
        slv.iterate();
        stressValues.push_back(slv.calculateStress(0));
        slv.writeVtkOutput(std::string("strain=") + std::to_string(cur_strain) + std::string("."));
        std::cout << std::endl;
    }

    std::string delim = "Stress values: ";
    for(auto val : stressValues) {
        std::cout << delim << val;
        delim = ", ";
    } std::cout << std::endl;

    /// write the final configuration to file
    std::cout << std::endl << "Deployment complete, writing output file" << std::endl;

    return EXIT_SUCCESS;
}

