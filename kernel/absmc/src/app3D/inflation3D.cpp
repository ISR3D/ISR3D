#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <chrono>

#include "absmc3D.h"

using namespace absmc;
using namespace absmc::graphics;

using util::createFileName;
using util::parseFileName;
using util::config;
using util::Timer;

using std::cout;
using std::cerr;
using std::endl;

typedef AgentBase<3> agent_t;
typedef Point<3, double> point_t;

/// main
int main(int argc, char *argv[])
{
    if (argc!=2) {
        cout << "Usage: deployStent configDir/configName" << endl;
        return EXIT_FAILURE;
    }

    /// read config filename; output to the config directory
    const std::string configFileName(argv[1]);
    const std::string configDir = parseFileName(configFileName).directory;
    const std::string outputDir = createFileName(configDir, "data.inflation", "");

    config().readFile(configFileName);

    /// agent input files
    const std::string inFileName = config().getValue<std::string>("deploy_input_file");
    const std::string baseName = parseFileName(inFileName).baseName;

    /// characteristic length scale
    const double L = config().getValue<double>("smc_mean_rad");

    /// load agents from file; includes vessel wall agents
    AgentContainer<3> agents;
    AgentFileReader<3>::readFile(createFileName(configDir, "stage2."+baseName, "dat"), AgentFactory3D(), agents);
    cout << agents.count() << " agents loaded." << endl;
    if (agents.count() == 0) {
        cout << "No vessel agents available, terminating." << endl;
        return 1;
    }

    /// iteration parameters
    const double iterTime = config().getValue<double>("deploy_iter_time");
    const int maxIter     = config().getValue<int>("deploy_max_iter");
    const int vtkIter     = config().getValue<int>("deploy_vtk_iter");
    const int datIter     = config().getValue<int>("deploy_dat_iter");

    /// initialise neighbourhoods
    VerletNeighbourDetector<agent_t> nbDetector(2 * L, 1.5*2*L); //works in general case, large neighbourhoods
    //nbDetector.init();

    /// Read bonds
    const std::string bondFileName = createFileName(configDir, "stage2."+baseName+"_nb", "dat");
    NBReader<3>::readFile(bondFileName, agents, true, true);

    /// set agents at longitudinal (x) domain boundaries immobile;
    const double lX = config().getValue<double>("lX");
    //const double lY = config().getValue<double>("lY");
    const double mX = config().getValue<double>("deploy_boundary_mobility_x");
    const double mY = config().getValue<double>("deploy_boundary_mobility_y");
    const double mZ = config().getValue<double>("deploy_boundary_mobility_z");
    agents.forAll(PSelectCloseToBoundary<agent_t>(0, 0.0, lX, L*4), FSetMobility<agent_t>(point_t(mX, mY, mZ) ) );
    size_t nAgents = agents.count();


    // SURFACE DETECTION
    auto t1 = std::chrono::system_clock::now();
    cout << endl << "Starting surface detection..." << endl;
    auto temp = std::chrono::system_clock::to_time_t(t1);
    cout << std::put_time(std::localtime(&temp), "%F %T") << endl;
    NeighbourCache<agent_t>* nbCacheForCiRule = new NeighbourCache<agent_t > (agents.getAgentVector());
    FDetectSurfaceFromNeighbours<agent_t> surfaceDetectorFunctor(L * 7, 200, &nbCacheForCiRule);
    agents.forAll(surfaceDetectorFunctor);
    auto t2 = std::chrono::system_clock::now();
    cout << "Surface detection complete, duration (s): " << std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << endl;

    size_t surfaceCount = 0;
    /// Count surface, set IEL and EEL
    for (size_t iAgent = 0; iAgent < nAgents; iAgent++) {
        agent_t* agent = agents.getAgentVector()[iAgent];
        if(agent->isSurface()){
            surfaceCount++;
        }
    }
    cout << surfaceCount << " surface agents" << endl;

    /// set rules for SMC agents (necrosis rule)
    const double smcMaxStress = config().getValue<double>("smc_max_stress");
    SMCNecrosisRule3D smcNecrosisRule(smcMaxStress);

    CompositeRule<SMC3D> smcRule;
    smcRule.add(&smcNecrosisRule);

    agents.forAll(tSMC3D, FSetAgentRule<SMC3D>(&smcRule) );

    /// set rules for IEL agents (iel breaking rule)
    const double ielMaxStrain = config().getValue<double>("iel_max_strain");
    const double ielMaxStress = config().getValue<double>("iel_max_stress");
    IELBreakingRule3D ielBreakingRule(ielMaxStrain , ielMaxStress);
    agents.forAll(tIEL3D, FSetAgentRule<IEL3D>(&ielBreakingRule) );

    const double yOffset = config().getValue<double>("deploy_y_offset");

//    /// stent deployment parameters
//    const double gap = config().getValue<double>("deploy_gap");
//    const double deploymentDepth = config().getValue<double>("deployment_depth");
//    const double dr = (deploymentDepth + gap) / maxIter;
//    cout << "Stent dr per macro iter = " << dr << endl;
//    double deltaR = - gap;

    ZeroUnaryForce3D zeroForce;
    ZeroBinaryForce3D zeroBinForce;
    // double unaryMagnitude = 7.794e-5; // 7.794e-5 = inflation test to 1 atm for r = 0.015 mm
    double unaryMagnitudeInflation = 1.039e-5; // 1.039e-5 = inflation test to 100 mmHg (13 333 Pa) for r = 0.015 mm
    // double unaryMagnitude = dr / iterTime; // stent displacement
    UnaryInflationForce displForce(unaryMagnitudeInflation, yOffset, 0.0);

    NeohookeanRepulsionForce3D bForce(L);
    Linear4thPowerAttractionForce3D smcAttrForce(L); /// Holzapfel, 0 strain

    LinearAttractionForce3D linAttrForce(0.001 * 5000); /// Averaged force for animal scenarios

    /// Forces for three-layer model
    Poly6thPowerAttractionForce3D adventitiaForce(0.001 * 1811.52 * 0.3536 * 2.0,
                                                  0.001 * 15262.42 * 0.3536 * 2.0,
                                                  0.001 * 0.34 * 0.3536 * 2.0,
                                                  0.001 * 9822726.48 * 0.3536 * 2.0,
                                                  0.001 * 125735191.07 * 0.3536 * 2.0,
                                                  0.001 * 38678553.78 * 0.3536 * 2.0);
    Poly6thPowerAttractionForce3D mediaForce(0.001 * 558 * 0.3536 * 2.0,
                                             0.,
                                             0.,
                                             0.,
                                             0.,
                                             0.001 * 91968100 * 0.3536 * 2.0);
    LinearAttractionForce3D neointimaForce(0.96);

    TypeSpecificBinaryForce threeLayerForce(&zeroBinForce);
    threeLayerForce.setDefaultType(tSMC3D); /// This type is used to calculate the characteristic interaction force
    threeLayerForce.add(&adventitiaForce, tEEL3D, tEEL3D); /// EEL agents are a stand-in for adventitia here, and IEL for intima
    threeLayerForce.add(&mediaForce, tSMC3D, tSMC3D);
    threeLayerForce.add(&neointimaForce, tIEL3D, tIEL3D);
    threeLayerForce.add(&mediaForce, tSMC3D, tEEL3D);
    threeLayerForce.add(&mediaForce, tSMC3D, tIEL3D);

    BinaryForce< AgentBase<3> > * binForce;
    binForce = &bForce;

    UnaryForce< AgentBase<3> > * uForce;
    uForce = &displForce;

    BinaryForce< AgentBase<3> > * bondForce;
    //bondForce = &smcAttrForce;
    bondForce = &linAttrForce;

    /// system-characteristic force
    const double compressionMagnitude = std::abs(binForce->calculateForce(L, L, L));
    const double stretchMagnitude = std::abs(bondForce->calculateForce(L, L, 2.3 * L));
    const double binaryMagnitude = std::max(compressionMagnitude, stretchMagnitude);
    const double charForce = std::max(binaryMagnitude, unaryMagnitudeInflation); //REAL UNITS

//    const double solverMaxDispl = (forceType == wmlc) ? 0.0005*L : 0.03*L; //REAL UNITS //0.001 for wmlc
    const double solverMaxDispl = 0.03*L; //0.03L is 1:20 per 30 iter //REAL UNITS //0.001 for wmlc
    const double maxDt = solverMaxDispl / charForce;
    const double solverMaxDt = 100. * maxDt;

    cout << "maxDt             = " << maxDt << endl;
    cout << "solverMaxDt       = " << solverMaxDt << endl;
    cout << "solverMaxDispl    = " << solverMaxDispl << endl;

    std::ofstream integratorLogFile(createFileName(outputDir, "fr", "dat").c_str() );

    RungeKuttaIntegratorMPI<agent_t>
            integrator(*uForce, *binForce, *bondForce, L, solverMaxDispl, solverMaxDt, nbDetector, integratorLogFile);

    FixedIntervalController fixedIntervalController(iterTime);

    /// desired convergence level
    const double eps      = config().getValue<double>("deploy_convergence_level");
    /// Convergence level for the final eqilibration
    const double epsFinal      = config().getValue<double>("deploy_convergence_level_final");

    cout << "charForce         = " << charForce << endl;
    cout << "eps               = " << eps << endl;
    ForceResidualMaxNormController frMaxNormController(std::numeric_limits<double>::infinity(), charForce, eps);
    ForceResidualMaxNormController frMaxNormControllerFinal(std::numeric_limits<double>::infinity(), charForce, epsFinal);

    std::cout << std::endl << "Interaction and bond forces for various distances" << std::endl;
    for (double dist = 1.0; dist < 3.0; dist += 0.1) {
        std::cout << "Distance = " << dist << " R; Interaction force = " << binForce->calculateForce(L, L, L*dist) <<
                                                "; Bond force = " << bondForce->calculateForce(L, L, L*dist) << std::endl;
    }

    /// vtp output parameters
    const std::string vtpScalars = config().getValue<std::string>("deploy_vtp_scalars");
    const std::string vtpVectors = config().getValue<std::string>("deploy_vtp_vectors");

    Timer tIter;
    Timer tTotal;


    /// steps of 25 mmhg
    for(size_t step = 1; step <= 4; step++){
        displForce.setMagnitude(unaryMagnitudeInflation * step / 4.);

        for (int iter=0; iter<maxIter; iter++) {
            if (vtkIter!=0 && iter%vtkIter==0) {

                tIter.stop();
                tTotal.stop();

                VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(),
                                                     createFileName(outputDir, baseName+"_"+std::to_string(step), "vtp", iter, 6),
                                                     vtpScalars, vtpVectors);

                cout << "iter=" << iter << ",  t(" << vtkIter << " it)="
                     << tIter.printElapsedTime() << ",  t(total)=" << tTotal.printElapsedTime()
                     << ", nSMC=" << agents.count(tSMC3D) << ", nIEL=" << agents.count(tIEL3D) << std::endl;

                tIter.start();
            }

            if (datIter!=0 && iter%datIter==0) {
                agents.writeFile(createFileName(outputDir, baseName+"_"+std::to_string(step), "dat", iter, 6) );
                NBWriter<3>::writeFile(agents, createFileName(outputDir, baseName+"_"+std::to_string(step)+"_nb", "dat", iter, 6) );
            }

            /// exec agent rules (iel breaking, smc necrosis)
            nbCacheForCiRule = new NeighbourCache<agent_t > (agents.getAgentVector());
            agents.execAgentRules();
            delete nbCacheForCiRule;

            /// main loop, smooth deployment: stent displacement is performed at each integration step
            integrator.integrate(agents.getAgentVector(), fixedIntervalController);

        }

        /// Final iteration with a smaller eps to bring the system closer to the equilibrium
        // uForce = &zeroForce;
        FixedIntervalController finFixedIntervalController(iterTime * 30);

        RungeKuttaIntegratorMPI<agent_t>
            finalIntegrator(*uForce, *binForce, *bondForce, L, solverMaxDispl / 2, solverMaxDt * 10, nbDetector, integratorLogFile);

        for (int iter=0; iter<5; iter++) {

            tIter.stop();
            tTotal.stop();

            VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(),
                                                 createFileName(outputDir, baseName, "vtp", iter, 6),
                                                 vtpScalars, vtpVectors);

            cout << "iter=" << iter << ",  t(" << vtkIter << " it)="
                 << tIter.printElapsedTime() << ",  t(total)=" << tTotal.printElapsedTime()
                 << ", nSMC=" << agents.count(tSMC3D) << ", nIEL=" << agents.count(tIEL3D) << std::endl;

            tIter.start();

            finalIntegrator.integrate(agents.getAgentVector(), finFixedIntervalController);
        }
    }

    /// write the final configuration to file
    cout << endl << "Deployment complete, writing output file" << endl;
    agents.writeFile(createFileName(configDir, "stage3."+baseName, "dat") );
    NBWriter<3>::writeFile(agents, createFileName(configDir, "stage3."+baseName+"_nb", "dat") );
    VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(), createFileName(configDir, "stage3."+baseName, "vtp"), vtpScalars, vtpVectors);
    cout << "Output successfully written" << endl;

    return EXIT_SUCCESS;
}
