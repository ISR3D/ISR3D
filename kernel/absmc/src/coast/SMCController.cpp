#include <stdexcept>
#include <random>
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
#include "absmc3D.h"

using libmuscle::Data;
using libmuscle::Instance;
using libmuscle::Message;
using ymmsl::Operator;

using namespace absmc;
using namespace absmc::graphics;

using util::logger;
using util::Random;
using util::createFileName;
using util::parseFileName;
using util::Timer;
//using muscle::cxa;

typedef AgentBase < 3 > agent_t;
typedef Point < 3, double> point_t;

/// Export cell list
/// Each data point is a list(list(x,y,z), r, typeId)
void exportCellListISR3D(AgentContainer < 3 > const& agentContainer, Data& data)
{
    const size_t nDim = 3;
    typedef AgentBase<nDim> agent_t;
    typedef Point<nDim, double> point_t;
    std::vector<agent_t*> const& agents = agentContainer.getAgentVector();

    size_t nAgents = agents.size();
    std::vector<double> posX(nAgents), posY(nAgents), posZ(nAgents), rad(nAgents);
    std::vector<int64_t> agentType(nAgents);

    for (std::size_t i = 0; i < agents.size(); ++i) {
        point_t pos = agents[i]->getPos();
        int64_t typeId = static_cast<int> (agents[i]->getTypeId());

        posX[i] = pos[0];
        posY[i] = pos[1];
        posZ[i] = pos[2];
        rad[i] = agents[i]->getR();
        agentType[i] = typeId;
//        data[i] = Data::list(
//            Data::list(pos[0], pos[1], pos[2]),
//            agents[i]->getR(),
//            typeId);
    }

    auto posX_result = Data::grid(posX.data(), {nAgents}, {"x"});
    auto posY_result = Data::grid(posY.data(), {nAgents}, {"y"});
    auto posZ_result = Data::grid(posZ.data(), {nAgents}, {"z"});
    auto rad_result = Data::grid(rad.data(), {nAgents}, {"r"});
    auto type_result = Data::grid(agentType.data(), {nAgents}, {"type"});

    data = Data::list(posX_result, posY_result, posZ_result, rad_result, type_result);

}

/// Import externally computed per-cell data.
void importExternalFields(AgentContainer < 3 > & agentContainer,
        const std::vector<double>& dataWssMax
        //vector<double> dataNeoDist
        )
{
    const size_t nDim = 3;
    typedef AgentBase<nDim> agent_t;

    std::vector<agent_t*> const& agents = agentContainer.getAgentVector();
    const size_t nAgents = agents.size();

    assert(dataWssMax.size() == nAgents);
//    assert(dataNeoDist.size() == nAgents);

    logger() << "Number of agents is: " << nAgents << std::endl;

    /// Read all values in sequence and assign them to agents
    for (size_t iAgent = 0; iAgent < nAgents; iAgent++)
    {
        if ((agents[iAgent]->getTypeId() == tEndo3D) ||
            (agents[iAgent]->getTypeId() == tSMC3D) ||
            (agents[iAgent]->getTypeId() == tIEL3D))
        {
            CellBase3D * agent;
            agent = (CellBase3D*) agents[iAgent];
            agent->setWssMax(dataWssMax[iAgent]);
            //agent->setDrugConc(dataDrugConc[iAgent]);
            // TODO: to be removed completely, unless we couple it to a different dd model
            agent->setDrugConc(0.);
            //agent->setNeoDist(dataNeoDist[iAgent]);
        }

        if (agents[iAgent]->getTypeId() == tObstacle3D)
        {
            Obstacle3D* agent = (Obstacle3D*) agents[iAgent];
            agent->setWssMax(dataWssMax[iAgent]);
            agent->setDrugConc(0.);
        }
    }
}

/// This function adds fenestrations to the IEL, replacing IEL agents with SMC agents with a set [probability]
void makeFenestrations (AgentContainer < 3 > &agentContainer, double probability) {
    const size_t nDim = 3;
    typedef AgentBase<nDim> agent_t;
    std::vector<agent_t*> & agents = agentContainer.getAgentVector();
    const size_t nAgents = agents.size();
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0., 1.);
    for (size_t iAgent = 0; iAgent < nAgents; iAgent++) {
        if (agents[iAgent]->getTypeId() == tIEL3D) {
            /// Change to SMC with a set probability
            if (distribution(generator) <= probability) {
                agentContainer.replace(iAgent, new SMC3D(agents[iAgent]->getPos(), agents[iAgent]->getPos0(), agents[iAgent]->getR()));
            }
        }
    }
}

/// Calculates reendothelialisation probability per hour for the moment [timeSinceStart] hours after stent implantation
double calculateEndothelialProbability(const double timeSinceStart, const double endPoint) {
    double increase_rate = (1-0.59)/24/(endPoint-3.00);         //59% endothelium after day 3 and 100% after 15 days

    /// Alternative, possibly more accurate reendothelization scenario. Little tested, results very similar to the currently used scenario.
    //double increase_rate = (1-0.95)/24/(endPoint-5.00);

    //double increase_rate = 1.0/24/endPoint;         //0%  endothelium after day 0 and 100% after 23 days. equation y=(100/23)x + 0. so we just take the slope, staring point is 0%.
    double endothelialProbability = 0.0;

    /// Alternative, possibly more accurate reendothelization scenario. Little tested, results very similar to the currently used scenario.
    /*if (timeSinceStart <= 120) {
        increase_rate = (0.95-0.53)/24/(endPoint-3.00);
    }*/

    if (timeSinceStart <= 72) {
        /// Alternative, possibly more accurate reendothelization scenario. Little tested, results very similar to the currently used scenario.
        //increase_rate = 0.53/72.0;

        increase_rate = 0.59/72.0; //going from 0 to 59% between days 0 an 3 (72 hours). Note: it could very well be that coverage does not end up to be exactly 59% after day 3.
    }

    //If the endpoint is not reached. EC probability will be 100% after "endPoint" days
    if (timeSinceStart <= 24*endPoint) {

        /// Alternative, possibly more accurate reendothelization scenario. Little tested, results very similar to the currently used scenario.
        /*if (timeSinceStart <= 120) {
            endothelialProbability = increase_rate/(1.0 - 0.53 - (timeSinceStart-72.0)*increase_rate);

            if (timeSinceStart <= 72) {
                endothelialProbability = increase_rate/(1.0 - timeSinceStart*increase_rate);
            }
        }
        else {
            endothelialProbability = increase_rate/(1.0 - 0.95 - (timeSinceStart-120.0)*increase_rate);
        }*/

        if (timeSinceStart <= 72) {
            endothelialProbability = increase_rate/(1.0 - timeSinceStart*increase_rate);
        } else {
            endothelialProbability = increase_rate/(1.0 - 0.59 - (timeSinceStart-72.0)*increase_rate);
        }
    } else {
        endothelialProbability = 1.0;
    }

    return endothelialProbability;
}


int main(int argc, char *argv[])
{
    try
    {
        Instance instance(argc, argv, {
                {Operator::S, {"cellMaxStressIn"}},
                {Operator::O_I, {"cellPositionsOut"}}
        });

        while (instance.reuse_instance()) {

            const std::string outPath = instance.get_setting_as<std::string> ("muscle_data_out_dir");
            const std::string logName = instance.get_setting_as<std::string> ("log_file");
            if(!logName.empty())
                logger().setFile(outPath + logName);
            logger() << "SMC native code started." << std::endl;
            logger() << "SMC writing output to " << outPath << std::endl;

            // agent input files
            const std::string inFileName = instance.get_setting_as<std::string> ("run_input_file");
            const std::string Stage3Path = instance.get_setting_as<std::string> ("muscle_data_in_dir");
            const std::string baseName = parseFileName(inFileName).baseName;

            const std::string centerlineFileName = instance.get_setting_as<std::string>("centerline_file");
            Centerline<3> centerline;
            centerline.readFromFile(createFileName(Stage3Path,centerlineFileName));

            // characteristic length scale
            const double L = instance.get_setting_as<double>("smc_mean_rad");

            // load agents from file; includes smc+iel+obstacle agents
            AgentContainer < 3 > agents;
            const std::string cellFileName = createFileName(Stage3Path, "stage3." + baseName, "dat");
            logger() << "Reading initial cells from file " << cellFileName << std::endl;
            AgentFileReader < 3 > ::readFile(cellFileName, AgentFactory3D(), agents);
            logger() << agents.count() << " agents loaded." << std::endl;
            logger() << "Total agents: " << agents.count(tSMC3D) << " SMCs, "
                                         << agents.count(tIEL3D) << " IEL agents, "
                                         << agents.count(tEEL3D) << " EEL agents." << std::endl;

            if(agents.count() == 0) {
                logger() << "NO AGENTS FOUND, ABORTING" << std::endl;
                return EXIT_FAILURE;
            }

            /// dummy WSS data for testing MUSCLE3
            std::vector<double> dummyWSS(agents.count(), 0.);
            auto dummyWSSData = Data::grid(dummyWSS.data(), {dummyWSS.size()}, {"id"});


            /// initialise neighbourhoods
    //        logger::finer("Initialising neighbourhoods and bonds...");
            VerletNeighbourDetector<agent_t> nbDetector(2 * L, 1.5*2*L); //works in general case, large neighbourhoods

            /// Read bonds
            const std::string bondFileName = createFileName(Stage3Path, "stage3."+baseName+"_nb", "dat");
            NBReader<3>::readFile(bondFileName, agents, true, true);
            logger() << "Neighbourhoods and bonds initialised." << std::endl;

            // Trim longitudinal boundaries to get nice flow inlet/outlet
            const double trimDistance = instance.get_setting_as<double>("boundary_trim_distance");
            agents.remove(PSelectCloseToBoundary<agent_t>(0, agents.getMinCoord(0), agents.getMaxCoord(0), trimDistance));

            logger() << "After trimming longitudinal boundaries: " << agents.count() << " agents left." << std::endl;
            logger() << "Total agents: " << agents.count(tSMC3D) << " SMCs, "
                                         << agents.count(tIEL3D) << " IEL agents, "
                                         << agents.count(tEEL3D) << " EEL agents." << std::endl;


            // set agents at longitudinal (x) domain boundaries immobile;
            const double BR = instance.get_setting_as<double>("Boundary_Inactivity_Range");
            const double mX = instance.get_setting_as<double>("run_boundary_mobility_x");
            const double mY = instance.get_setting_as<double>("run_boundary_mobility_y");
            const double mZ = instance.get_setting_as<double>("run_boundary_mobility_z");
            agents.forAll(PSelectCloseToBoundary<agent_t>(0, agents.getMinCoord(0), agents.getMaxCoord(0), BR),
                          FSetMobility<agent_t>(point_t(mX, mY, mZ) ) );

            // make fenestrations in the IEL
            logger() << "Making fenestrations..." << std::endl;
            const double fenProb = instance.get_setting_as<double>("fenestration_probability");
            makeFenestrations(agents, fenProb);

            /// read reendothelization endpoint
            const double endoEndPoint = instance.get_setting_as<double>("endo_endpoint");

            /// Set EEL immobile (surrounding tissue is stiff)
            agents.forAll(tEEL3D, FSetMobility<agent_t>(point_t(0.0, 0.0, 0.0) ) );

            /// set rules for SMC agents; no rules for IEL agents.
            logger() << "Setting SMC rules..." << std::endl;

            const double smcMaxStress = instance.get_setting_as<double>("smc_max_stress");
            SMCNecrosisRule3D smcNecrosisRule(smcMaxStress);

            const double ciWeightSMC = instance.get_setting_as<double>("ci_weight_smc");
            const double ciWeightIEL = instance.get_setting_as<double>("ci_weight_iel");
            const double ciWeightEEL = instance.get_setting_as<double>("ci_weight_eel");
            const double ciWeightObstacle = instance.get_setting_as<double>("ci_weight_obstacle");
            const double ciRangeFactor = instance.get_setting_as<double>("ci_range_factor");
            const double ciThresholdCount = instance.get_setting_as<double>("ci_threshold_count");
            NeighbourCache<agent_t>* nbCacheForCiRule = 0;
            SMCContactInhibitionRule3D smcContactInhibitionRule(ciWeightSMC, ciWeightIEL, ciWeightEEL, ciWeightObstacle, ciRangeFactor, ciThresholdCount, &nbCacheForCiRule);

            const double dt = instance.get_setting_as<double>("dt");
            double initialEndothelialProbability = 0.0;
            SMCNitricOxideRule3D smcNitricOxideRule3D(initialEndothelialProbability, dt);

            const double drugConcThreshold = instance.get_setting_as<double>("smc_drug_conc_threshold");
            SMCCellCycleRule3D smcCellCycleRule(drugConcThreshold, dt);

            ECMProductionRule3D ecmProductionRule3D(0.1, centerline);

    //        HMOX1Rule3D hmox1rule3D(dt);

            const double smcMaxStrain = instance.get_setting_as<double>("smc_max_strain");
            SMCBondBreakingRule3D smcBondBreakingRule(smcMaxStrain);

            CompositeRule<SMC3D> smcRule;
    //        smcRule.add(&smcNecrosisRule);
            smcRule.add(&smcContactInhibitionRule);
            smcRule.add(&smcNitricOxideRule3D);
            smcRule.add(&smcCellCycleRule);
            smcRule.add(&ecmProductionRule3D);
            smcRule.add(&smcBondBreakingRule);
    //        smcRule.add(&hmox1rule3D);
            agents.forAll(PSelectAgentsByType<agent_t> (tSMC3D), FSetAgentRule<SMC3D> (&smcRule));

            // set SMC agents at longitudinal (x) domain boundaries biologically inactive
            agents.forAll(PSelectCloseToBoundary<agent_t > (0, agents.getMinCoord(0), agents.getMaxCoord(0), BR), FSetSMC3DInactive());

            //set endothelium to denuded if the cell is under the balloon
            logger() << "Denuding endothelium under the balloon..." << std::endl;
            const double balloonExtension =  instance.get_setting_as<double>("balloon_extension");
            const double stentMin = agents.getMinCoord(0, tObstacle3D);
            const double stentMax = agents.getMaxCoord(0, tObstacle3D);
            logger() << "Stent from " << stentMin << " to " << stentMax << std::endl;
            agents.forAll(PSelectInterval<agent_t > (0, stentMin - balloonExtension, stentMax + balloonExtension), FDenudeSMC3DEndothelium());

            const double unaryMagnitude = 0.;
            ZeroUnaryForce3D zeroForce;
            ZeroBinaryForce3D zeroBForce;

            NeohookeanRepulsionForce3D bForce(L);
            SMCAttractionForce3D smcAttrForce(L);
            LinearAttractionForce3D linAttrForce(0.001 * 5000); /// Averaged force for animal scenarios

            BinaryForce< AgentBase<3> > * binForce;
            binForce = &bForce;

            UnaryForce< AgentBase<3> > * uForce;
            uForce = &zeroForce;

            BinaryForce< AgentBase<3> > * bondForce;
            bondForce = &linAttrForce;


            //CompositeUnaryForce uForce;
            //uForce.clear();

            /// system-characteristic force
            const double compressionMagnitude = std::abs(binForce->calculateForce(L, L, L));
            const double stretchMagnitude = std::abs(bondForce->calculateForce(L, L, 2.3 * L));
            const double binaryMagnitude = std::max(compressionMagnitude, stretchMagnitude);
            const double charForce = std::max(binaryMagnitude, unaryMagnitude); //REAL UNITS
            logger() << "charForce         = " << charForce << std::endl;

            const double solverMaxDispl = 0.03*L; //REAL UNITS
            /// desired convergence level
            const double eps = instance.get_setting_as<double>("run_convergence_level");

            /// maximal time step (from rough stability estimation)
            const double maxDt = solverMaxDispl / charForce;
            logger() << "MaxDtReal        = " << std::endl;
            /// maximal time step passed to solver
            const double solverMaxDt = 50 * maxDt;

            logger() << "dt [hours]        = " << dt << std::endl;
            logger() << "maxDt             = " << maxDt << std::endl;
            logger() << "solverMaxDt       = " << solverMaxDt << std::endl;
            logger() << "solverMaxDispl    = " << solverMaxDispl << std::endl;
            logger() << "charForce         = " << charForce << std::endl;
            logger() << "eps               = " << eps << std::endl;

            logger() << std::endl << "Interaction and bond forces for various distances" << std::endl;
            for (double dist = 1.0; dist < 3.0; dist += 0.1) {
                logger() << "Distance = " << dist << " R; Interaction force = " << binForce->calculateForce(L, L, L*dist) <<
                                                        "; Bond force = " << bondForce->calculateForce(L, L, L*dist) << std::endl;
            }

            std::ofstream integratorLogFile(createFileName(outPath, "fr", "dat").c_str());
            RungeKuttaIntegratorMPI<agent_t>
                integrator(*uForce, *binForce, *bondForce, L, solverMaxDispl, solverMaxDt, nbDetector, integratorLogFile);

            ForceResidualMaxNormController frMaxNormController(std::numeric_limits<double>::infinity(), charForce, eps);
            //ForceResidualMaxNormController frMaxNormController(solverMaxDt * 500, charForce, eps); //arbitrary limit of steps

            /// iteration parameters
            const int vtkIter = instance.get_setting_as<int64_t>("run_vtk_iter");
            const int datIter = instance.get_setting_as<int64_t>("run_dat_iter");

            /// vtp output parameters
            const std::string vtpScalars = instance.get_setting_as<std::string>("run_vtp_scalars");
            const std::string vtpVectors = instance.get_setting_as<std::string>("run_vtp_vectors");

            Timer tIter;
            Timer tTotal;
            int iter = 0;

            /// main loop
            logger() << "SMC starting main loop." << std::endl;

            double newEndothelialProbability = 0.0;
            Random::init(0); /// trying to control the randomness by enforcing repetition

            const double timeMax = instance.get_setting_as<double>("max_time");
            double timeInHours = 0.; //time since the simulation started

            while (timeInHours < timeMax)
            {

                tIter.start();

                Data cellsOut;
                exportCellListISR3D(agents, cellsOut);
                instance.send("cellPositionsOut", Message(timeInHours, cellsOut));

                logger() << "Starting integrator" << std::endl;
                integrator.integrate(agents.getAgentVector(), frMaxNormController);
                logger() << "Integrator exited" << std::endl;

                tIter.stop();
                tTotal.stop();
                logger() << "Iteration compute walltime: " << tIter.printElapsedTime() << ", total compute walltime: " << tTotal.printElapsedTime() << std::endl;

                logger() << "Simulated time since agent-based simulation start = " << timeInHours << " hours, before receiving data from conduits." << std::endl;

                //Message dummyWSSMessage(timeInHours, dummyWSSData);
                auto wssIn = instance.receive("cellMaxStressIn");//, dummyWSSMessage);
                auto wss_data_ptr = wssIn.data().elements<double>();
                std::vector<double> dataWssMax (wss_data_ptr, wss_data_ptr + wssIn.data().size());

                logger() << "received WSS" << std::endl;
    //            vector<double> dataNeoDist = muscle::env::receiveDoubleVector("CellNeointimaDistance");
    //            logger::info("received NeoDist");
                importExternalFields(agents, dataWssMax);//, dataNeoDist);

                if ((vtkIter != 0) && (iter % vtkIter == 0)) {
                    logger() << "Writing VTK file" << std::endl;
                    VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(),
                                                         createFileName(outPath, baseName, "vtp", iter, 6),
                                                         vtpScalars, vtpVectors);

                    logger() << "iter=" << iter << ", nSMC=" << agents.count(tSMC3D) << ", nIEL=" << agents.count(tIEL3D)
                             << ", nEEL=" << agents.count(tEEL3D) << ", nObstacles=" << agents.count(tObstacle3D) << std::endl;
                }

                if (datIter!=0 && iter%datIter==1) {
                    logger() << "Writing checkpoint file" << std::endl;
                    agents.writeFile(createFileName(outPath, "stage4."+baseName, "dat", iter, 4) );
                    NBWriter<3>::writeFile(agents, createFileName(outPath, "stage4."+baseName+"_nb", "dat", iter, 4) );
                }

                logger() << "SMC: timeInHours=" << timeInHours << ",  before executing biology rules" << std::endl;
                // set SMC agents at longitudinal (x) domain boundaries biologically inactive
                agents.forAll(PSelectCloseToBoundary<agent_t > (0, agents.getMinCoord(0), agents.getMaxCoord(0), BR), FSetSMC3DInactive());

                // exec agent rules
                nbCacheForCiRule = new NeighbourCache<agent_t > (agents.getAgentVector());

                newEndothelialProbability = calculateEndothelialProbability(iter * dt, endoEndPoint);
                smcNitricOxideRule3D.setNewEndothelialProbability(newEndothelialProbability);
                if (iter > 0)
                    agents.execAgentRules();

                delete nbCacheForCiRule;

                iter++;
                timeInHours += dt;
            }

            /// write the final configuration to file
            agents.writeFile(createFileName(outPath, baseName, "dat", iter, 6));
            VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(), createFileName(outPath, "stage4." + baseName, "vtp"), vtpScalars, vtpVectors);

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
