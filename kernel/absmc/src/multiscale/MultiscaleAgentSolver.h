#ifndef ABSMC_MULTISCALE_AGENT_SOLVER_H
#define ABSMC_MULTISCALE_AGENT_SOLVER_H

#include <libmuscle/libmuscle.hpp>
#include <ymmsl/ymmsl.hpp>
#include <random>
#include <iterator>
#include <algorithm>
#include "core/agentSolver.h"
#include "absmc3D.h"

using libmuscle::Data;
using libmuscle::Instance;
using libmuscle::Message;
using ymmsl::Operator;

namespace absmc {

template<size_t nDim>
class MultiscaleAgentSolver : public AgentSolver<nDim> {

public:

    using typename AgentSolver<nDim>::agent_t;
    using typename AgentSolver<nDim>::point_t;
    using typename AgentSolver<nDim>::metrics_t;

    MultiscaleAgentSolver(Instance& instance_) :
    instance(&instance_) {

        outputDir = instance->get_setting_as<std::string> ("muscle_data_out_dir");
        logName = instance->get_setting_as<std::string> ("log_file");
        if(!logName.empty())
            logger().setFile(outputDir + logName);
        logger() << "SMC native code started." << std::endl;
        logger() << "SMC writing output to " << outputDir << std::endl;

        timeMax = instance->get_setting_as<double>("max_time");
        dt = instance->get_setting_as<double>("dt");
        logger() << "dt [hours]        = " << dt << std::endl;

        /// read reendothelization endpoint
        endoEndPoint = instance->get_setting_as<double>("endo_endpoint");
    }

    /// Export cell list
    /// Each data point is a list(list(x,y,z), r, typeId)
    void exportCellListISR3D(Data& data)
    {
        std::vector<agent_t*> const& agentsVector = agents.getAgentVector();

        nAgents = agents.count();
        std::vector<double> posX(nAgents), posY(nAgents), posZ(nAgents), rad(nAgents);
        std::vector<int64_t> agentType(nAgents);

        for (std::size_t i = 0; i < agentsVector.size(); ++i) {
            point_t pos = agentsVector[i]->getPos();
            int64_t typeId = static_cast<int> (agentsVector[i]->getTypeId());

            posX[i] = pos[0];
            posY[i] = pos[1];
            posZ[i] = pos[2];
            rad[i] = agentsVector[i]->getR();
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
    void importExternalFields(const std::vector<double>& dataWssMax
            //vector<double> dataNeoDist
            )
    {
        std::vector<agent_t*> const& agentsVector = agents.getAgentVector();
        nAgents = agentsVector.size();

        assert(dataWssMax.size() == nAgents);
    //    assert(dataNeoDist.size() == nAgents);

        logger() << "Number of agents is: " << nAgents << std::endl;

        /// Read all values in sequence and assign them to agents
        for (size_t iAgent = 0; iAgent < nAgents; iAgent++)
        {
            if ((agentsVector[iAgent]->getTypeId() == tEndo3D) ||
                (agentsVector[iAgent]->getTypeId() == tSMC3D) ||
                (agentsVector[iAgent]->getTypeId() == tIEL3D) ||
                (agentsVector[iAgent]->getTypeId() == tECM3D))
            {
                CellBase3D * agent;
                agent = (CellBase3D*) agentsVector[iAgent];
                agent->setWssMax(dataWssMax[iAgent]);
                //agent->setDrugConc(dataDrugConc[iAgent]);
                // TODO: to be removed completely, unless we couple it to a different dd model
                agent->setDrugConc(0.);
                //agent->setNeoDist(dataNeoDist[iAgent]);
            }

            if (agentsVector[iAgent]->getTypeId() == tObstacle3D)
            {
                Obstacle3D* agent = (Obstacle3D*) agentsVector[iAgent];
                agent->setWssMax(dataWssMax[iAgent]);
                agent->setDrugConc(0.);
            }
        }
    }


    void setStage4Rules() {
        /// set rules for SMC agents; no rules for IEL agents.
        logger() << "Setting SMC rules..." << std::endl;

        const double smcMaxStress = instance->get_setting_as<double>("smc_max_stress");
        smcNecrosisRule = std::make_shared< SMCNecrosisRule3D >(smcMaxStress);

        const double ciWeightSMC = instance->get_setting_as<double>("ci_weight_smc");
        const double ciWeightIEL = instance->get_setting_as<double>("ci_weight_iel");
        const double ciWeightEEL = instance->get_setting_as<double>("ci_weight_eel");
        const double ciWeightObstacle = instance->get_setting_as<double>("ci_weight_obstacle");
        const double ciRangeFactor = instance->get_setting_as<double>("ci_range_factor");
        const double ciThresholdCount = instance->get_setting_as<double>("ci_threshold_count");
        smcContactInhibitionRule = std::make_shared<SMCContactInhibitionRule3D> (ciWeightSMC, ciWeightIEL, ciWeightEEL, ciWeightObstacle, ciRangeFactor, ciThresholdCount, &nbCacheForCiRule);

        smcNitricOxideRule3D = std::make_shared<SMCNitricOxideRule3D> ();

        const double drugConcThreshold = instance->get_setting_as<double>("smc_drug_conc_threshold");
        smcCellCycleRule = std::make_shared<SMCCellCycleRule3D> (drugConcThreshold, dt);

        const double ecmLooseThreshold = instance->get_setting_as<double>("ecm_loose_threshold");
        const double ecmProductionProbability = instance->get_setting_as<double>("ecm_prod_probability");
        ecmProductionRule3D = std::make_shared<ECMProductionRule3D> (ecmLooseThreshold, ecmProductionProbability, centerline);

        //HMOX1Rule3D hmox1rule3D(dt);

        const double smcMaxStrain = instance->get_setting_as<double>("smc_max_strain");
        smcBondBreakingRule = std::make_shared<SMCBondBreakingRule3D> (smcMaxStrain);

        smcRule = std::make_shared< CompositeRule<SMC3D> > ();
//        smcRule->add(smcNecrosisRule);
        smcRule->add(smcContactInhibitionRule);
        smcRule->add(smcNitricOxideRule3D);
        smcRule->add(smcCellCycleRule);
        smcRule->add(ecmProductionRule3D);
        smcRule->add(smcBondBreakingRule);
//        smcRule.add(&hmox1rule3D);

        agents.forAll(tSMC3D, FSetAgentRule<SMC3D>(smcRule) );
    }

    /// Sets surface SMC's as EC-covered to obtain at least "precentage" coverage
    /// This version uses WSS value to determine surface cells, alternative: CI rule or direct surface detection
    /// The cells that are set as inactive are excluded from the calcultion
    void setSurfaceECCoveredSMCPercentage(double percentage) {
        // count(SMC, WSS != 0., BActive)

        std::vector<SMC3D*> covered;
        std::vector<SMC3D*> uncovered;

        std::vector<agent_t*> const& agentsVector = agents.getAgentVector();

        for(auto& agentPtr : agentsVector) {
            if(agentPtr->getTypeId() == tSMC3D) {
                SMC3D* smcPtr = dynamic_cast<SMC3D*>(agentPtr);
                if(smcPtr->getBaFlag() && (smcPtr->getWssMax() != 0.)) {
                    if(smcPtr->getSelectedFlag()) {
                        covered.push_back(smcPtr);
                    } else {
                        uncovered.push_back(smcPtr);
                    }
                }
            }
        }

        size_t surfaceTotal = covered.size() + uncovered.size();
        int sampleSize = std::max((int)std::round(surfaceTotal * percentage - covered.size()), 0);
        std::vector<SMC3D*> toCover;

        logger() << "On the vessel surface, there are " << uncovered.size() << " SMCs not covered by endothelium (" <<
                 (double)(uncovered.size()) / (double) (surfaceTotal + 1.) * 100.  << "% of the wound area)." << std::endl;

        if(sampleSize == 0) return;

        /// draw from the list until percentage is met
        std::sample(uncovered.begin(), uncovered.end(), std::back_inserter(toCover),
                sampleSize, std::mt19937{std::random_device{}()});
        for(auto& agentPtr : toCover) {
            agentPtr->setSelectedFlag(true);
        }
    }

    // iterates in dt increments until timeMax is reached
    void iterateOverTime() {
        Timer tIter;
        Timer tTotal;
        /// main loop
        logger() << "SMC starting main loop." << std::endl;

        int iter = 0;
        Random::init(0); /// trying to control the randomness by enforcing repetition
        double timeInHours = 0.; //time since the simulation started

        while (timeInHours < timeMax)
        {
            tIter.start();

            Data cellsOut;
            exportCellListISR3D(cellsOut);
            instance->send("cellPositionsOut", Message(timeInHours, cellsOut));

            logger() << "Starting integrator" << std::endl;
            integrator->integrate(agents.getAgentVector(), *controller);
            logger() << "Integrator exited" << std::endl;

            tIter.stop();
            tTotal.stop();
            logger() << "Iteration compute walltime: " << tIter.printElapsedTime() << ", total compute walltime: " << tTotal.printElapsedTime() << std::endl;

            logger() << "Simulated time since agent-based simulation start = " << timeInHours << " hours, before receiving data from conduits." << std::endl;

            //Message dummyWSSMessage(timeInHours, dummyWSSData);
            auto wssIn = instance->receive("cellMaxStressIn");//, dummyWSSMessage);
            auto wss_data_ptr = wssIn.data().elements<double>();
            std::vector<double> dataWssMax (wss_data_ptr, wss_data_ptr + wssIn.data().size());

            logger() << "received WSS" << std::endl;
            importExternalFields(dataWssMax);

            if (vtkIter!=0 && iter%vtkIter==0) {
                tIter.stop();
                tTotal.stop();
                logger() << "Writing VTK file" << std::endl;
                graphics::VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(),
                                                               createFileName(outputDir, baseName, "vtp", iter, 6),
                                                               vtpScalars, vtpVectors);
                logger() << "iter=" << iter << ",  t(" << vtkIter << " it)="
                          << tIter.printElapsedTime() << ",  t(total)=" << tTotal.printElapsedTime()
                          << ", nSMC=" << agents.count(tSMC3D) << ", nIEL=" << agents.count(tIEL3D)
                          << ", nEEL=" << agents.count(tEEL3D) << ", nObstacles=" << agents.count(tObstacle3D) << std::endl;
                tIter.start();
            }
            if (datIter!=0 && iter%datIter==0) {
                logger() << "Writing checkpoint file" << std::endl;
                agents.writeFile(createFileName(outputDir, "stage4."+baseName, "dat", iter, 6) );
                NBWriter<nDim>::writeFile(agents, createFileName(outputDir, "stage4."+baseName+"_nb", "dat", iter, 6));
            }

            logger() << "SMC: timeInHours=" << timeInHours << ",  before executing biology rules" << std::endl;
            // set SMC agents at longitudinal (x) domain boundaries biologically inactive
            this->setSMCsInactive(instance->get_setting_as<double>("balloon_axial_protrusion"));

            // update endothelium coverage
            setSurfaceECCoveredSMCPercentage(this->calculateEndothelialPercentage(timeInHours, endoEndPoint));

            // exec agent rules
            nbCacheForCiRule = new NeighbourCache<agent_t> (agents.getAgentVector());
            if (iter > 0)
                agents.execAgentRules();
            delete nbCacheForCiRule;

            iter++;
            timeInHours += dt;
        }
    }

private:
    using AgentSolver<nDim>::agents;
    using AgentSolver<nDim>::nAgents;

    using AgentSolver<nDim>::centerline;

    using AgentSolver<nDim>::solverName;
    using AgentSolver<nDim>::configDir;
    using AgentSolver<nDim>::dataInputDir;
    using AgentSolver<nDim>::outputDir;

    using AgentSolver<nDim>::inFileName;
    using AgentSolver<nDim>::outFileName;
    using AgentSolver<nDim>::baseName;
    using AgentSolver<nDim>::vtpScalars;
    using AgentSolver<nDim>::vtpVectors;
    using AgentSolver<nDim>::L;
    using AgentSolver<nDim>::sigma;

    using AgentSolver<nDim>::binForce;
    using AgentSolver<nDim>::bondForce;
    using AgentSolver<nDim>::uForce;
    using AgentSolver<nDim>::integrator;
    using AgentSolver<nDim>::controller;

    using AgentSolver<nDim>::maxIter;
    using AgentSolver<nDim>::vtkIter;
    using AgentSolver<nDim>::datIter;

    using AgentSolver<nDim>::charForce;
    using AgentSolver<nDim>::binaryMagnitude;
    using AgentSolver<nDim>::unaryMagnitude;
    using AgentSolver<nDim>::solverMaxDispl;
    using AgentSolver<nDim>::solverMaxDt;

    using AgentSolver<nDim>::smcNecrosisRule;
    using AgentSolver<nDim>::smcContactInhibitionRule;
    using AgentSolver<nDim>::smcNitricOxideRule3D;
    using AgentSolver<nDim>::smcCellCycleRule;
    using AgentSolver<nDim>::ecmProductionRule3D;
    using AgentSolver<nDim>::smcBondBreakingRule;

    using AgentSolver<nDim>::smcRule;
    using AgentSolver<nDim>::ielBreakingRule;


    Instance* instance;
    double dt;
    double timeMax;
    double endoEndPoint;

    std::string logName;

    NeighbourCache<agent_t>* nbCacheForCiRule;

};

}

#endif
