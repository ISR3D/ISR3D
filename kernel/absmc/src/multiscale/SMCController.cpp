#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <limits>
#include <math.h>
#include <libmuscle/libmuscle.hpp>
#include <ymmsl/ymmsl.hpp>
#include <stdlib.h>
#include "multiscale/MultiscaleAgentSolver.h"

using libmuscle::Data;
using libmuscle::Instance;
using libmuscle::Message;
using ymmsl::Operator;

using namespace absmc;
using namespace absmc::graphics;

int main(int argc, char *argv[])
{
    try
    {
        Instance instance(argc, argv, {
                {Operator::S, {"cellMaxStressIn"}},
                {Operator::O_I, {"cellPositionsOut"}}
        });

        MultiscaleAgentSolver<3> slv(instance);

        while (instance.reuse_instance()) {
            slv.setDataInputDir(instance.get_setting_as<std::string> ("muscle_data_in_dir"));
            slv.setInFileName(instance.get_setting_as<std::string> ("run_input_file"));
            slv.setBaseNameFromIn();

            // characteristic length scale
            slv.setCharLengthScales(instance.get_setting_as<double>("smc_mean_rad"), 0.);

            slv.readDat3DFile("stage3." + slv.getBaseName() + ".dat");
            slv.initVerletNeighbourDetector();
            slv.readBondsFromFile("stage3." + slv.getBaseName() + "_nb.dat");
            logger() << "Neighbourhoods and bonds initialised." << std::endl;

            slv.readCenterlineFromFile(instance.get_setting_as<std::string>("centerline_file"));

            /// dummy WSS data for testing MUSCLE3
//            std::vector<double> dummyWSS(agents.count(), 0.);
//            auto dummyWSSData = Data::grid(dummyWSS.data(), {dummyWSS.size()}, {"id"});

            // Trim longitudinal boundaries to get nice flow inlet/outlet
            slv.trimBoundaries(instance.get_setting_as<double>("boundary_trim_distance"));

            slv.setBoundaryMobility(instance.get_setting_as<double>("run_boundary_mobility_x"),
                                    instance.get_setting_as<double>("run_boundary_mobility_y"),
                                    instance.get_setting_as<double>("run_boundary_mobility_z"),
                                    instance.get_setting_as<double>("Boundary_Inactivity_Range"));

            /// Set EEL immobile (surrounding tissue is stiff)
            slv.setEelImmobile();

            // make fenestrations in the IEL
            logger() << "Making fenestrations..." << std::endl;
            slv.makeFenestrations(instance.get_setting_as<double>("fenestration_probability"));

            /// set rules for SMC agents; no rules for IEL agents.
            slv.setStage4Rules();
            slv.setBalloonDenudation(instance.get_setting_as<double>("balloon_axial_protrusion"));
            slv.setSMCsInactive(instance.get_setting_as<double>("balloon_axial_protrusion"));

            slv.setForces(&slv.bForce, &slv.linAttrForce, &slv.zeroForce);
            slv.calculateCharForce();

            slv.initializeIntegrator();

            /// desired convergence level
            slv.setForceResidualMaxNormController(instance.get_setting_as<double>("run_convergence_level"));

            /// iteration parameters
            slv.setVtkIter(instance.get_setting_as<int64_t>("run_vtk_iter"));
            slv.setDatIter(instance.get_setting_as<int64_t>("run_dat_iter"));

            /// vtp output parameters
            slv.setVtpParameters(instance.get_setting_as<std::string>("run_vtp_scalars"),
                                 instance.get_setting_as<std::string>("run_vtp_vectors"));

            slv.iterateOverTime();

            logger() << "SMC model finished iterating, writing final output." << std::endl;

            /// write the final configuration to file
            slv.writeDatOutputNoBonds("stage4.");
            slv.writeVtkOutput("stage4.");

            //////////////////////////////////////////////////
            logger() << "SMC finishing native code." << std::endl;

        }
        return EXIT_SUCCESS;
    } catch (std::runtime_error& e)
    {
        logger() << "RUNTIME ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...)
    {
        logger() << "UNKNOWN ERROR" << std::endl;
        return EXIT_FAILURE;
    }
}
