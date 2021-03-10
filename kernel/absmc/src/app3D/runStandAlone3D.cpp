#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <limits>
#include <boost/timer.hpp>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <omp.h>

#include "absmc3D.h"

using namespace absmc;
using namespace absmc::graphics;

using util::createFileName;
using util::parseFileName;
using util::config;

using std::cout;
using std::cerr;
using std::endl;

typedef AgentBase<3> agent_t;
typedef Point<3, double> point_t;

class FSetState : public std::binary_function<size_t, agent_t, void> {
public:
    void operator()(size_t id, agent_t* agentBase) {
        SMC3D* agent = dynamic_cast<SMC3D*>(agentBase);
        if (agent) { agent->setState(SMC3D::SG2M); }
    }
};

/// This executable is outdated. Use at your own risk
/// main
int main(int argc, char *argv[])
{
    if (argc!=3) {
        cout << "Usage: runStandAlone configDir outputRoot" << endl;
        return EXIT_FAILURE;
    }

    const std::string configDir(argv[1]);
    const std::string outputRoot(argv[2]);
    if (mkdir(outputRoot.c_str(),0755)==-1 && errno != EEXIST) {
    cout<<"miserably failed to create outputRoot <" << outputRoot << ">: "<<((std::string)strerror(errno)) << endl;
    return 1;
    }

    if (mkdir((outputRoot+"/data.runStandAlone").c_str(),0755)==-1 && errno != EEXIST) {cout<<"miserably failed to create data.runStandAlone <" << outputRoot << "/data.runStandAlone>: "<<((std::string)strerror(errno)) << endl; return 1; }
    const std::string outputDir = createFileName(outputRoot, "data.runStandAlone", "");

    const std::string configFileName  = createFileName(configDir, "absmc", "cfg");
    config().readFile(configFileName);

    // agent input files
    const std::string inFileName = config().getValue<std::string>("run_input_file");
    const std::string baseName = parseFileName(inFileName).baseName;

    // morse potential parameters (dimensionless)
    const double U1 = config().getValue<double>("U1");
    const double r1 = config().getValue<double>("r1");
    const double U2 = config().getValue<double>("U2");
    const double r2 = config().getValue<double>("r2");
    // equilibrium distance (dimensionless)
    const double r0 = r1*r2/(r2-r1)*log(U1/U2*r2/r1);
    // characteristic length scale
    const double L = config().getValue<double>("smc_mean_rad");
    // boundary range
    const double BR = config().getValue<double>("Boundary_Inactivity_Range");

    // load agents from file; includes smc+iel+obstacle agents
    AgentContainer<3> agents;
    AgentFileReader<3>::readFile(createFileName(configDir, "stage3."+baseName, "dat"), AgentFactory3D(), agents);
    cout << agents.count() << " agents loaded." << endl;

    // initialise neighbourhoods
    AccumulativeNeighbourDetector<agent_t> nbDetector(r0*L);
    const std::string nbFileName = createFileName(configDir, "stage3."+baseName+"_nb", "dat");
    NBReader<3>::readFile(nbFileName, agents);

    // set agents at longitudinal (x) domain boundaries immobile;
    const double lX = config().getValue<double>("lX");
    const double mX = config().getValue<double>("run_boundary_mobility_x");
    const double mY = config().getValue<double>("run_boundary_mobility_y");
    const double mZ = config().getValue<double>("run_boundary_mobility_z");
    agents.forAll(PSelectCloseToBoundary<agent_t>(0, 0.0, lX, BR), FSetMobility<agent_t>(point_t(mX, mY, mZ) ) );

    //ADDED set mobility of EEL AND SMCS to zero
    agents.forAll(tEEL3D, FSetMobility<agent_t>(point_t(0.0, 0.0, 0.0) ) );
    agents.forAll(tIEL3D, FSetMobility<agent_t>(point_t(0.0, 0.0, 0.0) ) );
//    agents.forAll(tSMC3D, FSetMobility<agent_t>(point_t(0.0, 0.0, 0.0) ) );

    // set SMC agents at longitudinal (x) domain boundaries and at the outer surface biologically inactive
    agents.forAll(PSelectCloseToBoundary<agent_t>(0, 0.0, lX, BR), FSetSMC3DInactive() );
    // set rules for SMC agents; no rules for IEL agents.
    const double smcMaxStress = config().getValue<double>("smc_max_stress");
    SMCNecrosisRule3D smcNecrosisRule(smcMaxStress);

    const double ciWeightSMC          = config().getValue<double>("ci_weight_smc");
    const double ciWeightIEL          = config().getValue<double>("ci_weight_iel");
    const double ciWeightEEL          = config().getValue<double>("ci_weight_eel");
    const double ciWeightObstacle     = config().getValue<double>("ci_weight_obstacle");
    const double ciRangeFactor        = config().getValue<double>("ci_range_factor");
    const double ciThresholdCount     = config().getValue<double>("ci_threshold_count");
    const double ciRange = ciRangeFactor*r0*L;
    NeighbourCache<agent_t>* nbCacheForCiRule = 0;
    SMCContactInhibitionRule3D smcContactInhibitionRule(ciWeightSMC, ciWeightIEL, ciWeightEEL, ciWeightObstacle, ciRange, ciThresholdCount, &nbCacheForCiRule);

    const double drugConcThreshold    = config().getValue<double>("smc_drug_conc_threshold");
    const double wssMaxThreshold      = config().getValue<double>("smc_wss_max_threshold");
    SMCCellCycleRule3D smcCellCycleRule(drugConcThreshold, wssMaxThreshold);

    //pretending we have a blood flow code running


    //    logger::info("Sizes nagents=%d, OSI=%d, WSS=%d, Drug=%d",nAgents,dataWssOsi.size(),dataWssMax.size(),dataDrugConc.size());
        int timeStep = 0;
    //new Nitric oxide rule
    SMCNitricOxideRule3D smcNitricOxideRule3D(0.0); //is timestep defined?

    CompositeRule<SMC3D> smcRule;
//    smcRule.add(&smcNecrosisRule);
    smcRule.add(&smcContactInhibitionRule);
    smcRule.add(&smcCellCycleRule);
    smcRule.add(&smcNitricOxideRule3D);

    // set rules for bulk smc agents only; other smc agents are left biologically inactive with nil agent rule.
    const double outerRadius = config().getValue<double>("outer_radius");
    agents.forAll(PSelectBulkAgentsOfSpecifiedType<agent_t>(0.0, lX, 0.0, 0.0, outerRadius, 4.0*0.5*r0*L, tSMC3D),
                  FSetAgentRule<SMC3D>(&smcRule) );

    //set initial number of SMCS in NO rule
    smcNitricOxideRule3D.setInitial_Num_SMCs(agents.count(tSMC3D));
/*
    //////////////// physics parameters for Morse force (dimensionless)
    // cell-cell stiffness at equilibrium distance is zero, so kEq is at 0.001
    //    const double young = 0.1; //MPa
    const double youngeq = 0.77;//0.5*young/(1.0-0.27*0.27);
    const double B = 1.0/8.3439;//20.0*youngeq/(9.0*math::pi); //MPa = N/mm2
    const double kEq = 0.00273753;//0.5*B*math::pi*0.000716223/L;// N/mm/L
    // maximal time step (from rough stability estimation)
    const double maxDt = 2.0 / kEq; // L*mm/N
    // maximal time step passed to solver
    const double solverMaxDt = maxDt;  // mm/N
    // maximal displacement
    const double solverMaxDispl = 0.05*2.0*L;// mm //0.1*0.5*r0;
*/

    //////////////// physics parameters for Morse force (dimensionless)
    // cell-cell stiffness at equilibrium distance
    const double kEq = 8.0;//U1/(r1*r1)*exp(-r0/r1) - U2/(r2*r2)*exp(-r0/r2);
    // maximal time step (from rough stability estimation)
    const double maxDt = 2.0 / kEq;;
    // maximal time step passed to solver
    const double solverMaxDt = 0.01*maxDt;;
    // maximal displacement
    const double solverMaxDispl = 0.01*0.5*r0;

/*
    //////////////// physics parameters for Morse force (dimensionless)
    // cell-cell stiffness at equilibrium distance
    const double B = 1.0/8.3439;
    const double kEq = B*0.0027/0.015/0.015;//0.01/6.0/0.015;//U1/(r1*r1)*exp(-r0/r1) - U2/(r2*r2)*exp(-r0/r2);
    // maximal time step (from rough stability estimation)
    const double maxDt = 2.0 / kEq;;
    // maximal time step passed to solver
    const double solverMaxDt = 0.1*maxDt;;
    // maximal displacement
    const double solverMaxDispl = 0.1*0.5*r0;
*/

    cout << "U1                = " << U1 << endl;
    cout << "r1                = " << r1 << endl;
    cout << "U2                = " << U2 << endl;
    cout << "r2                = " << r2 << endl;
    cout << "r0                = " << r0 << endl;
    cout << "kEq               = " << kEq << endl;;
    cout << "maxDt             = " << maxDt << endl;
    cout << "solverMaxDt       = " << solverMaxDt << endl;
    cout << "solverMaxDispl    = " << solverMaxDispl << endl;

    ZeroUnaryForce3D uForce;
    Bilinear3D bForce(U1, r1, U2, r2, L);

    std::ofstream integratorLogFile(createFileName(outputDir, "fr", "dat").c_str() );
    AdaptiveEulerIntegrator<agent_t, ZeroUnaryForce3D, Bilinear3D>
            integrator(uForce, bForce, L, solverMaxDispl, solverMaxDt, nbDetector, integratorLogFile);

    // iteration parameters
    const double iterTime = config().getValue<double>("run_iter_time");
    const int maxIter     = config().getValue<int>("run_max_iter");
    const int vtkIter     = config().getValue<int>("run_vtk_iter");
    const int datIter     = config().getValue<int>("run_dat_iter");
    FixedIntervalController fixedIntervalController(iterTime);

    // system-characteristic force
    const double charForce = kEq * r0 * 0.1;
    // desired convergence level
    const double eps      = config().getValue<double>("run_convergence_level");

    cout << "charForce         = " << charForce << endl;
    cout << "eps               = " << eps << endl;
    //cahrForce*eps = 0.000108268 for morse seems like 0.00005 for beyondhertz (already multiplied by 0.015, still need to apply B).
    ForceResidualMaxNormController frMaxNormController(std::numeric_limits<double>::infinity(), charForce, eps);

    // vtp output parameters
    const std::string vtpScalars = config().getValue<std::string>("run_vtp_scalars");
    const std::string vtpVectors = config().getValue<std::string>("run_vtp_vectors");

    boost::timer timer;
    double tIter = 0.0;
    double tTotal = 0.0;
    double omptimefinal = 0.0;

    double omptime = omp_get_wtime();

    // main loop

    timeStep = 0;
    int newSMCs = 0;
    smcNitricOxideRule3D.setNewSMCCounter(newSMCs);

    int oldNumOfSMCs = agents.count(tSMC3D);
    int newNumOfSMCs = oldNumOfSMCs;
    int selectedSMCs = 0;
    int initial_num_SMCs = agents.count(tSMC3D);
    double new_endothelial_probability = 0.0;
    util::Random::init(0); //trying to control the randomness by allowing repetition


    for (int iter=0; iter<maxIter; iter++) {

        const size_t nAgents = agents.getAgentVector().size();
        //pretending we have a blood flow code running
        for (size_t iAgent = 0; iAgent < nAgents; iAgent++)
        {
            if ((agents.getAgentVector()[iAgent]->getTypeId() == tSMC3D) || (agents.getAgentVector()[iAgent]->getTypeId() == tIEL3D))
            {
                CellBase3D * agent;
                agent = (CellBase3D*) agents.getAgentVector()[iAgent];
                //        agent->setWssOsi(dataWssOsi[iAgent]);
                agent->setWssMax(12.0);

            }
        }


        if (iter%vtkIter==0) {

            tIter = timer.elapsed();
            tTotal += tIter;
            omptimefinal = - omptime + omp_get_wtime();

            VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(),
                                                 createFileName(outputDir, baseName, "vtp", iter, 6),
                                                 vtpScalars, vtpVectors);

            cout << "iter=" << iter << ",  t(" << vtkIter << " it)=";
            cout << std::setprecision(4) << tIter << ",  t(total)=" << tTotal;
            cout << ",  t(omp)=" << std::setprecision(4) << omptimefinal;
            cout << ", nSMC=" << agents.count(tSMC3D) << ", nIEL=" << agents.count(tIEL3D)<< ", nEEL=" << agents.count(tEEL3D)<< ", nSEL=" << selectedSMCs
            << ", nNEW=" << newSMCs << ", prob=" << 100.0*new_endothelial_probability << std::endl;

            timer.restart();
        }

        if (datIter!=0 && iter%datIter==0) {
            NBWriter<3>::writeFile(agents, createFileName(configDir, "stage4."+baseName+"_nb", "dat") );
            agents.writeFile(createFileName(outputDir, baseName, "dat", iter, 6) );
        }
         // set SMC agents at longitudinal (x) domain boundaries and at the outer surface biologically inactive
        agents.forAll(PSelectCloseToBoundary<agent_t>(0, 0.0, lX, BR), FSetSMC3DInactive() );
        // exec agent rules
        nbCacheForCiRule = new NeighbourCache<agent_t>(agents.getAgentVector() );



        double increase_rate = 0.0205/24.0;                 // considering 59% endothelium after day 3rd and 100% after 23 days. equation y=(41/20)x + 52.85. so we just take the slope, staring point is 59%.
        //       double increase_rate = 0.41/12.0/24.0;             // considering 59% endothelium after day 3rd and 100% after 15 days. equation y=(41/12)x + 48.75. so we just take the slope, staring point is 59%.
        //double increase_rate = 1.0/23.0/24.0;         //0.181159420289855    // considering 0% endothelium after day 0 and 100% after 23 days. equation y=(100/23)x + 0. so we just take the slope, staring point is 0%.
        double old_endothelial_probability = 0.0;            // probability assuming a fixed number of cells. This is just used to calculate the probabilty defined in the next line.
        //endothelial_probability = 0.0;
        new_endothelial_probability = 0.0;            // probability that includes the newly produced SMCs, More the SMC divide into daugter cells, more the probabilty is.

        if (iter <=72) {
            increase_rate = 0.59/3.0/24.0; //going from 0 to 59% between days 0 an 3. Note: it could very well be that coverage does not end up to be exactly 59% after day 3.
        } else {
            increase_rate = 0.41/12.0/24.0;
        }

        if (iter <= 360) {                     // EC probabilty will be 100% after 15 days, 15x24=360
            //old_endothelial_probability = 100.0*increase_rate/(47.15 - iter*increase_rate); //NOW iter SHOULD BE IN HOURS, VALID FOR 3 to 23, we use 100-52.85=47.15
            //        old_endothelial_probability = 100.0*increase_rate/(51.25 - iter*increase_rate); //NOW iter SHOULD BE IN HOURS, VALID FOR 3 to 15, we use 100-48.75=51.25
            if (iter <=72) {
                old_endothelial_probability = increase_rate/(1.0 - iter*increase_rate); //0.001934235976789 //NOW iter SHOULD BE IN HOURS, VALID FOR 0 to 23
            } else {
                old_endothelial_probability = increase_rate/(1.0 - 72.0*0.59/3.0/24.0 - (iter-72.0)*increase_rate);
            }
            //double percent_cells_left = 47.15 - iter*increase_rate;        //VALID FOR 3 to 23, we use 100-52.85=47.15
            //        double percent_cells_left = 51.25 - iter*increase_rate;        //VALID FOR 3 to 15, we use 100-48.75=51.25
            double percent_cells_left = 1.0 - iter*increase_rate;    //0.936594202898551    //VALID FOR 0 to 23
            int num_cells_left = (percent_cells_left * initial_num_SMCs); //274991.550724637679798
            //        int newSMCs = 0;
            //    endothelial_probability = 100.0*(old_endothelial_probability*num_cells_left/100.0 + newSMCs)/(num_cells_left + newSMCs);
            new_endothelial_probability = old_endothelial_probability;//(old_endothelial_probability*num_cells_left + newSMCs)/(num_cells_left + newSMCs);  // this is required as the total number of SMCs are also increasing with time.

            //                if (iter == 72) { /*endothelial_probability = 59.0;*/new_endothelial_probability = 59.0;}     // first step for Nakazawa et al.
        } else { /*endothelial_probability = 100.0;*/new_endothelial_probability = 1.0; }


        //    printf("endprob = %f \n", new_endothelial_probability);


        if (new_endothelial_probability > 1.0) { /*endothelial_probability =100.0;*/new_endothelial_probability = 1.0; }


        smcNitricOxideRule3D.setnew_endothelial_probability(100.0*new_endothelial_probability);

        agents.execAgentRules();
        newNumOfSMCs = agents.count(tSMC3D);
        selectedSMCs = agents.countSelected(tSMC3D);
        newSMCs = newNumOfSMCs - oldNumOfSMCs;
        smcNitricOxideRule3D.setNewSMCCounter(newSMCs);
        delete nbCacheForCiRule;

        // ... and restore equilibrium
        nbDetector.updateSlowTimeScale(agents.getAgentVector() );
        integrator.integrate(agents.getAgentVector(), frMaxNormController);
    }

    // write the final configuration to file
    agents.writeFile(createFileName(configDir, "stage4."+baseName, "dat") );
     NBWriter<3>::writeFile(agents, createFileName(configDir, "stage4."+baseName+"_nb", "dat") );
    VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(), createFileName(configDir, "stage4."+baseName, "vtp"), vtpScalars, vtpVectors);

    return EXIT_SUCCESS;
}
