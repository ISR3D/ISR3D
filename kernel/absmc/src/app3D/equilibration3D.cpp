#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

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

/// Shuffle agent positions; the direction of the displacement vector is uniformly distributed,
/// and its norm is Gauss-distributed with specified mean and sigma. Only the specified fraction of agents is shuffled,
/// the remainder is left at their initial position and set immobile. This allows to fix a certain number of dof when
/// studying if the ensemble returns to a pre-shuffling equilibrium configuration.
/// Agents already marked as immobile are not displaced either; note that if their number is non-zero,
/// the fraction of agents actually displaced will differ from the specified fraction.
template<class Agent>
void shuffleAgents(std::vector<Agent*> & agents, double meanShift, double sigmaShift, double shuffleFraction)
{
    typedef Point<Agent::nDim, double> point_t;

    const size_t nAgents = agents.size();
    size_t nFixed = 0;
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {

        const double shift = util::Random::getNextGauss(meanShift, sigmaShift);
        const double theta = util::Random::getNext(math::pi);
        const double phi   = util::Random::getNext(2.0*math::pi);

        point_t newPos = agents[iAgent]->getPos();

        newPos[0] += shift * sin(theta) * cos(phi);
        newPos[1] += shift * sin(theta) * sin(phi);
        newPos[2] += shift * cos(theta);

        // displace only the specified fraction of agents, fix the remainder.
        // do not displace agents which were already fixed.
        const double p = util::Random::getNext(1.0);
        if (p < shuffleFraction && agents[iAgent]->getMobility() == point_t(1.0, 1.0, 1.0) ) {
            agents[iAgent]->setPos(newPos);
        }
        else {
            agents[iAgent]->setMobility(point_t(0.0, 0.0, 0.0) );
            nFixed++;
        }
    }
    std::cout << (nAgents-nFixed) << " agents displaced randomly with meanShift=" << meanShift << ", sigmaShift=" << sigmaShift << ". ";
    std::cout << nFixed << " agents undisplaced and set immobile." << std::endl;
}

typedef AgentBase<3> agent_t;
typedef Point<3, double> point_t;

/// main
int main(int argc, char *argv[])
{
    if (argc!=2) {
        cout << "Usage: equilibrate configDir/configName" << endl;
        return EXIT_FAILURE;
    }

    const std::string configFileName(argv[1]);
    const std::string configDir = parseFileName(configFileName).directory;
    const std::string outputDir = createFileName(configDir, "data.equilibration", "");

    config().readFile(configFileName);

    // agent input file
    const std::string inFileName = config().getValue<std::string>("eq_input_file");
    const std::string baseName = parseFileName(inFileName).baseName;

      // characteristic length scale
    const double L = config().getValue<double>("smc_mean_rad");
    const double sigma = config().getValue<double>("smc_rad_sigma");

    // load agents from file includes smc+iel+obstacle agents
    AgentContainer<3> agents;
    AgentFileReader<3>::readFile(createFileName(configDir, "stage1."+baseName, "dat"), AgentFactory3D(), agents);
    cout << agents.count() << " agents loaded." << endl;

    // create neighbour detector but do not initialise neighbourhoods...
    AccumulativeNeighbourDetector<agent_t> nbDetector(2.2 * L + 2 * sigma);
    // StaticNeighbourDetector<agent_t> nbDetector(2 * L + 2 * sigma);
    // ... instead load neighbourhoods from file and print nb statistics;
    const std::string nbFileName = createFileName(configDir, "stage1."+baseName+"_nb", "dat");
    NBReader<3>::readFile(nbFileName, agents);

    // set agents at longitudinal (x) domain boundaries immobile.
    const double mX = config().getValue<double>("eq_boundary_mobility_x");
    const double mY = config().getValue<double>("eq_boundary_mobility_y");
    const double mZ = config().getValue<double>("eq_boundary_mobility_z");
    agents.forAll(PSelectCloseToBoundary<agent_t>(0, agents.getMinCoord(0), agents.getMaxCoord(0), L),
                  FSetMobility<agent_t>(point_t(mX, mY, mZ) ) );


    // randomly displace agents; only the specified fraction will be displaced,
    // the remainder is left in position and set immobile
    const double meanShift  = config().getValue<double>("eq_mean_shift");
    const double sigmaShift = config().getValue<double>("eq_sigma_shift");
    const double shuffleFraction = config().getValue<double>("eq_shuffle_fraction");
    shuffleAgents(agents.getAgentVector(), meanShift, sigmaShift, shuffleFraction);

    ZeroUnaryForce3D zeroForce;
    ZeroBinaryForce3D zeroBinForce;

    NeohookeanRepulsionForce3D bForce(L);
    Linear4thPowerAttractionForce3D smcAttrForce(L); /// Holzapfel, 0 strain

    LinearAttractionForce3D linAttrForce(0.001 * 5000); /// Averaged force for animal scenarios

    /// Forces for three-layer model

    ///Adventitia, Cilla et al.
    Poly6thPowerAttractionForce3D adventitiaForce(0.07946 * 0.3536 * 2.0,
                                             0.,
                                             0.,
                                             0.,
                                             0.,
                                             1255.50532 * 0.3536 * 2.0);
    ///Media, Cilla et al.
    Poly6thPowerAttractionForce3D mediaForce(0.13897 * 0.3536 * 2.0,
                                             0.,
                                             0.,
                                             0.,
                                             0.,
                                             5.53748937e+02 * 0.3536 * 2.0);
    ///Intima, Chiastra et al.
    Poly6thPowerAttractionForce3D neointimaForce(0.0949584 * 0.3536 * 2.0,
                                             0.,
                                             0.,
                                             0.,
                                             80.9860722 * 0.3536 * 2.0,
                                             2.13579344 * 0.3536 * 2.0);

    TypeSpecificBinaryForce threeLayerForce(&zeroBinForce);
    threeLayerForce.setDefaultType(tSMC3D); /// This type is used to calculate the characteristic interaction force
    threeLayerForce.add(&adventitiaForce, tEEL3D, tEEL3D); /// EEL agents are a stand-in for adventitia here, and IEL for intima
    threeLayerForce.add(&mediaForce, tSMC3D, tSMC3D);
    threeLayerForce.add(&neointimaForce, tIEL3D, tIEL3D);
    threeLayerForce.add(&mediaForce, tSMC3D, tEEL3D);
    threeLayerForce.add(&mediaForce, tSMC3D, tIEL3D);


    BinaryForce< AgentBase<3> > * binForce;
    binForce = &bForce;

    BinaryForce< AgentBase<3> > * bondForce;
    bondForce = &linAttrForce;

    /// system-characteristic force
    const double compressionMagnitude = std::abs(binForce->calculateForce(L, L, L));
    const double stretchMagnitude = std::abs(bondForce->calculateForce(L, L, 2.3 * L));
    const double binaryMagnitude = std::max(compressionMagnitude, stretchMagnitude);
    //const double charForce = std::max(binaryMagnitude, unaryMagnitude); //REAL UNITS
    const double charForce = binaryMagnitude; //REAL UNITS

    // maximal displacement
    const double solverMaxDispl = 0.03*L; //REAL UNITS
    // maximal time step (from rough stability estimation)
    const double maxDt = solverMaxDispl / charForce;
    // maximal time step passed to solver
    const double solverMaxDt = 10. * maxDt;

    cout << "maxDt             = " << maxDt << endl;
    cout << "solverMaxDt       = " << solverMaxDt << endl;
    cout << "solverMaxDispl    = " << solverMaxDispl << endl;

    ZeroUnaryForce3D uForce;
    std::ofstream integratorLogFile(createFileName(outputDir, "fr", "dat").c_str() );

    RungeKuttaIntegratorMPI<agent_t>
            integrator(uForce, *binForce, *bondForce, L, solverMaxDispl, solverMaxDt, nbDetector, integratorLogFile);

    // iteration parameters
    const int maxIter     = config().getValue<int>("eq_max_iter");
    const int vtkIter     = config().getValue<int>("eq_vtk_iter");
    const int datIter     = config().getValue<int>("eq_dat_iter");

    // vtp output parameters
    const std::string vtpScalars = config().getValue<std::string>("eq_vtp_scalars");
    const std::string vtpVectors = config().getValue<std::string>("eq_vtp_vectors");

    std::ofstream l2DiffStr(createFileName(outputDir, "l2diff", "dat").c_str() );

    Timer tIter;
    Timer tTotal;

    // desired convergence level
    const double eps      = config().getValue<double>("eq_convergence_level");

    cout << "charForce         = " << charForce << endl;
    cout << "eps               = " << eps << endl;
    ForceResidualMaxNormController frMaxNormController(1000, charForce, eps);

    // main loop
    for (int iter=0; iter<maxIter; iter++) {

        if (vtkIter!=0 && iter%vtkIter==0) {

            tIter.stop();
            tTotal.stop();

            VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(),
                                                 createFileName(outputDir, baseName, "vtp", iter, 6),
                                                 vtpScalars, vtpVectors);

            double nbDistMean, nbDistStdDev;
            computeNbDistStat<agent_t>(agents, &nbDistMean, &nbDistStdDev);

            l2DiffStr << iter << "  \t"
                      << integrator.getTotalIter() << " \t"
                      << integrator.getTotalTime() << " \t"
                      << std::setprecision(10) << l2Difference(agents.getAgentVector()) << " \t"
                      << nbDistMean << " \t" << nbDistStdDev << endl;

            cout << "iter=" << iter << ",  t(" << vtkIter << " it)="
                 << tIter.printElapsedTime() << ",  t(total)=" << tTotal.printElapsedTime() << std::endl;

            tIter.start();
        }

        if (datIter!=0 && iter%datIter==0) {
            agents.writeFile(createFileName(outputDir, baseName, "dat", iter, 6) );
        }

        // create bonds from the detected neighbourhoods
        /// a bit dirty, can be a container-wide thing or pre-set in csv
        size_t nAgents = agents.count();
        for (size_t iAgent = 0; iAgent < nAgents; iAgent++) {
            agents.getAgentVector()[iAgent]->addNeighboursToBonds(true, 12); // max bonds set to 12 since that's the number of neighbours in a closest packing
        }
        cout << "Bonds initialized from neighbours" << endl;

        nbDetector.updateSlowTimeScale(agents.getAgentVector() );
        integrator.integrate(agents.getAgentVector(), frMaxNormController);
    }

    // set initial positions of cells to current positions
    agents.forAll(FResetPos0<agent_t>() );

    // write the final configuration to file
    agents.writeFile(createFileName(configDir, "stage2."+baseName, "dat") ); // agents
    NBWriter<3>::writeFile(agents, createFileName(configDir, "stage2."+baseName+"_nb", "dat"), true); // bonds
    VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(), createFileName(configDir, "stage2."+baseName, "vtp"), vtpScalars, vtpVectors); // vtk image

    return EXIT_SUCCESS;
}
