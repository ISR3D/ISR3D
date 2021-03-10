#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>

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

enum ForceType { hertz, wmlc };
/// main
int main(int argc, char *argv[])
{
    if (argc!=2) {
        cout << "Usage: deployStent configDir/configName" << endl;
        return EXIT_FAILURE;
    }

    std::cout << "Deployment simulation started at " << Timer::printTime() << std::endl;

    /// read config filename; output to the config directory
    const std::string configFileName(argv[1]);
    const std::string configDir = parseFileName(configFileName).directory;
    const std::string outputDir = createFileName(configDir, "data.deployment", "");

    config().readFile(configFileName);

    /// agent input files
    const std::string inFileName = config().getValue<std::string>("deploy_input_file");
    const std::string baseName = parseFileName(inFileName).baseName;
    const std::string stentFileName = config().getValue<std::string>("deploy_stent_file");
    const std::string balloonFileName = config().getValue<std::string>("deploy_balloon_file");

    /// characteristic length scale
    const double L = config().getValue<double>("smc_mean_rad");

    /// iteration parameters
    const double iterTime = config().getValue<double>("deploy_iter_time");
    const int maxIter     = config().getValue<int>("deploy_max_iter");
    const int vtkIter     = config().getValue<int>("deploy_vtk_iter");
    const int datIter     = config().getValue<int>("deploy_dat_iter");

    /// load agents from file; includes vessel wall agents. Obstacles are loaded by the stent class, see below.
    AgentContainer<3> agents;
    AgentFileReader<3>::readFile(createFileName(configDir, "stage2."+baseName, "dat"), AgentFactory3D(), agents);
    cout << agents.count() << " agents loaded." << endl;
    if (agents.count() == 0) {
        cout << "No vessel agents available, terminating." << endl;
        return 1;
    }

    /// initialise neighbourhoods
    VerletNeighbourDetector<agent_t> nbDetector(2 * L, 1.5*2*L); //works in general case, large neighbourhoods
    //nbDetector.init();

    /// Read bonds
    const std::string bondFileName = createFileName(configDir, "stage2."+baseName+"_nb", "dat");
    NBReader<3>::readFile(bondFileName, agents, true, true);

    /// set agents at longitudinal (x) domain boundaries immobile;
    const double lX = config().getValue<double>("lX");
    const double mX = config().getValue<double>("deploy_boundary_mobility_x");
    const double mY = config().getValue<double>("deploy_boundary_mobility_y");
    const double mZ = config().getValue<double>("deploy_boundary_mobility_z");
    agents.forAll(PSelectCloseToBoundary<agent_t>(0, 0.0, lX, L*4), FSetMobility<agent_t>(point_t(mX, mY, mZ) ) );

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
    size_t nAgents = agents.count();

    /// load stent from file
    const double yOffset = config().getValue<double>("deploy_y_offset");
    SimpleStent3D stent(createFileName(configDir, stentFileName), yOffset, 0.0 );
    stent.exportObstacles(agents);
    cout << stent.count() << " obstacles loaded, Y offset is " << yOffset << endl;

    /// stent deployment parameters
    const double gap = config().getValue<double>("deploy_gap");
    const double deploymentDepth = config().getValue<double>("deployment_depth");
    const double dr = (deploymentDepth + gap) / maxIter;
    cout << "Stent dr per macro iter = " << dr << endl;
    double deltaR = - gap;

    /// move stent to the center of the vessel and expand it to touch the inner radius of the vessel
    /// careful, there are already cell centers at inner_radius
    /// Can cause all kinds of problems for curved vessels. The correct procedure would be to detect the agents and to expand the stent to touch them.
    /// Might even be a good idea to use voxelization for collision detection.
    /// Set the last value to "true" if you intend to use this with post-deployment stent shapes, e.g. from MRI
    /// stent.affix(config().getValue<double>("inner_radius") - L - gap, config().getValue<double>("lX"), false);
    stent.moveToCenterX(config().getValue<double>("lX"));

    /// Generate simple balloon
    // const double balloonXMin = config().getValue<double>("deploy_balloon_x_min");
    // const double balloonXMax = config().getValue<double>("deploy_balloon_x_max");
    // cout << "Generating balloon between " << balloonXMin << " and "<< balloonXMax << endl;
    // SimpleBalloon3D balloon(agents, yOffset, 0.0, balloonXMin, balloonXMax, L * 2);

    /// load balloon from file
    SimpleBalloon3D balloon(createFileName(configDir, balloonFileName), yOffset, 0.0 );
    //balloon.scale(stent.getScalingValue(), config().getValue<double>("lX"));
    balloon.moveToCenterX(config().getValue<double>("lX"));
    balloon.exportObstacles(agents);

    cout << "Balloon and stent scaling factor is " << stent.getScalingValue() << endl;
    cout << balloon.count() << " balloon obstacles loaded, Y offset is " << yOffset << endl;
    cout << agents.count() << " agents in total" << endl;

    const std::string centerlineFileName = config().getValue<std::string>("centerline_file");
    Centerline<3> centerline;
    centerline.readFromFile(createFileName(configDir,centerlineFileName));


    ZeroUnaryForce3D zeroForce;
    double unaryMagnitude = dr / iterTime; // stent displacement
    UnaryStentDisplacementForce displForce(unaryMagnitude, centerline);

    NeohookeanRepulsionForce3D bForce(L);
    SMCAttractionForce3D smcAttrForce(L);
    LinearAttractionForce3D linAttrForce(0.001 * 5000); /// Averaged force for animal scenarios

    BinaryForce< AgentBase<3> > * binForce;
    binForce = &bForce;

    UnaryForce< AgentBase<3> > * uForce;
    uForce = &displForce;

    BinaryForce< AgentBase<3> > * bondForce;
    bondForce = &linAttrForce;

    /// system-characteristic force
    const double compressionMagnitude = std::abs(binForce->calculateForce(L, L, L));
    const double stretchMagnitude = std::abs(bondForce->calculateForce(L, L, 2.3 * L));
    const double binaryMagnitude = std::max(compressionMagnitude, stretchMagnitude);
    const double charForce = std::max(binaryMagnitude, unaryMagnitude); //REAL UNITS

//    const double solverMaxDispl = (forceType == wmlc) ? 0.0005*L : 0.03*L; //REAL UNITS //0.001 for wmlc
    const double solverMaxDispl = 0.05*L; //0.03L is 1:20 per 30 iter //REAL UNITS //0.001 for wmlc
    const double maxDt = solverMaxDispl / charForce;
    const double solverMaxDt = 500. * maxDt;

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

    cout << endl << "Interaction forces for various distances" << endl;
    for (double dist = 1.0; dist < 2.3; dist += 0.1) {
        cout << "Distance = " << dist << " R; Force = " << binForce->calculateForce(L, L, L*dist) << endl;
    }

    cout << endl << "Bond forces for various distances" << endl;
    for (double dist = 1.9; dist < 3.0; dist += 0.1) {
        cout << "Distance = " << dist << " R; Force = " << bondForce->calculateForce(L, L, L*dist) << endl;
    }

    /// vtp output parameters
    const std::string vtpScalars = config().getValue<std::string>("deploy_vtp_scalars");
    const std::string vtpVectors = config().getValue<std::string>("deploy_vtp_vectors");


//    /// SURFACE DEBUG
//    nbCacheForCiRule = new NeighbourCache<agent_t > (agents.getAgentVector());
//    DetectSurfaceFromNeighboursRule3D surfaceDetectorRule(L * 7, 200, &nbCacheForCiRule);
//    agents.forAll(FSetAgentRule<CellBase3D>(&surfaceDetectorRule) );
//    agents.execAgentRules();

//    // vtp output parameters
//    const std::string vtpScalarsTest = config().getValue<std::string>("reader_vtp_scalars");
//    const std::string vtpVectorsTest = config().getValue<std::string>("reader_vtp_vectors");

//    // write the final configuration to file
//    VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(), createFileName(configDir, "testSurfaceEq."+baseName, "vtp"), vtpScalarsTest, vtpVectorsTest);

    Timer tIter;
    Timer tTotal;

    for (int iter=0; iter<maxIter; iter++) {
        if (vtkIter!=0 && iter%vtkIter==0) {

            tIter.stop();
            tTotal.stop();

            VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(),
                                                 createFileName(outputDir, baseName, "vtp", iter, 6),
                                                 vtpScalars, vtpVectors);

            cout << "iter=" << iter << ",  t(" << vtkIter << " it)="
                 << tIter.printElapsedTime() << ",  t(total)=" << tTotal.printElapsedTime() << ", deltaR=" << deltaR
                 << ", nSMC=" << agents.count(tSMC3D) << ", nIEL=" << agents.count(tIEL3D) << std::endl;

            tIter.start();
        }

        if (datIter!=0 && iter%datIter==0) {
            agents.writeFile(createFileName(outputDir, baseName, "dat", iter, 6) );
            NBWriter<3>::writeFile(agents, createFileName(outputDir, baseName+"_nb", "dat", iter, 6) );
        }

        /// exec agent rules (iel breaking, smc necrosis)
        agents.execAgentRules();

        /// main loop, smooth deployment: stent displacement is performed at each integration step
        integrator.integrate(agents.getAgentVector(), fixedIntervalController);
        deltaR += dr;

    }

    cout << "Removing balloon" << endl;
    balloon.removeBalloon(agents);
    cout << "Balloon removed, equilibrating..." << endl;
    /// Final iteration with a smaller eps to bring the system closer to the equilibrium
    uForce = &zeroForce;
    FixedIntervalController finFixedIntervalController(iterTime * 20);

    RungeKuttaIntegratorMPI<agent_t>
        finalIntegrator(*uForce, *binForce, *bondForce, L, solverMaxDispl / 2, solverMaxDt * 10, nbDetector, integratorLogFile);

    for (int iter=0; iter<5; iter++) {
        tIter.stop();
        tTotal.stop();

        VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(),
                                             createFileName(outputDir, baseName, "vtp", iter, 6),
                                             vtpScalars, vtpVectors);

        cout << "iter=" << iter << ",  t(" << vtkIter << " it)="
             << tIter.printElapsedTime() << ",  t(total)=" << tTotal.printElapsedTime() << ", deltaR=" << deltaR
             << ", nSMC=" << agents.count(tSMC3D) << ", nIEL=" << agents.count(tIEL3D) << std::endl;

        tIter.start();

        finalIntegrator.integrate(agents.getAgentVector(), finFixedIntervalController);
    }

    /// write the final configuration to file
    cout << endl << "Deployment complete, writing output file" << endl;
    agents.writeFile(createFileName(configDir, "stage3."+baseName, "dat") );
    NBWriter<3>::writeFile(agents, createFileName(configDir, "stage3."+baseName+"_nb", "dat") );
    VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(), createFileName(configDir, "stage3."+baseName, "vtp"), vtpScalars, vtpVectors);
    cout << "Output successfully written" << endl;

    return EXIT_SUCCESS;
}

