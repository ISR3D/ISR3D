#ifndef ABSMC_AGENT_SOLVER_H
#define ABSMC_AGENT_SOLVER_H

// This class defines an object with common parameters and ,methods of AB simulations.
// E.g. parameter reading, interaction with the config, force instances and so on are defined here.
// Warning: using functions out of order (e.g. reading from a config before specifying its filename)
// is for the most part not handled and will crash your simulation.
// Possible future improvement: make child classes for specific simulations, read parameters in the constructor

#include <iostream>
#include <memory>
#include <random>
#include <util/zzInclude.h>

#include "core/agentBase.h"
#include "core/agentContainer.h"
#include "core/agentFactory.h"
#include "core/agentFunctors.h"
#include "core/agentPredicates.h"
#include "core/neighbourIO.h"
#include "core/neighbourDetector.h"
#include "core/staticNeighbourDetector.h"
#include "core/currentNeighbourDetector.h"
#include "core/verletNeighbourDetector.h"
#include "core/accumulativeNeighbourDetector.h"

#include "core/integrator.h"
#include "core/adaptiveRKIntegrator.h"
#include "core/integratorController.h"

#include "graphics/vtkWriter.h"
#include "model3D/zzInclude.h"

using util::createFileName;
using util::parseFileName;
using util::config;
using util::Timer;
using util::logger;
using util::Random;


namespace absmc {

template<size_t nDim>
class AgentSolver {

public:
    typedef AgentBase<nDim> agent_t;
    typedef Point<nDim, double> point_t;

    typedef EuclideanMetrics<nDim, double> metrics_t;


    // Specific force values used in simulations

    ZeroUnaryForce3D zeroForce;
    ZeroBinaryForce3D zeroBinForce;
    NeohookeanRepulsionForce3D bForce;
    //static const Linear4thPowerAttractionForce3D smcAttrForce(L); /// Holzapfel, 0 strain

    LinearAttractionForce3D linAttrForce; /// Averaged force for animal scenarios

    /// Forces for the three-layer model

    ///Adventitia, Cilla et al.
    Poly6thPowerAttractionForce3D adventitiaForce;
    ///Media, Cilla et al.
    Poly6thPowerAttractionForce3D mediaForce;
    ///Intima, Chiastra et al.
    Poly6thPowerAttractionForce3D neointimaForce;

    TypeSpecificBondForce threeLayerForce;

    // call this before using the three layer force
    void setupThreeLayerForce() {
        threeLayerForce.setDefaultType(tSMC3D); /// This type is used to calculate the characteristic interaction force
        threeLayerForce.add(&adventitiaForce, tEEL3D, tEEL3D); /// EEL agents are a stand-in for adventitia here, and IEL for intima
        threeLayerForce.add(&mediaForce, tSMC3D, tSMC3D);
        threeLayerForce.add(&neointimaForce, tIEL3D, tIEL3D);
        threeLayerForce.add(&mediaForce, tSMC3D, tEEL3D);
        threeLayerForce.add(&mediaForce, tSMC3D, tIEL3D);
    }

    AgentSolver(std::string const& solverName_="ABSMC solver") :
         linAttrForce(2.18),
         adventitiaForce(0.014,
                         0.,
                         0.,
                         0.,
                         40.,
                         175.),
         mediaForce(0.05,
                     0.,
                     0.,
                     0.,
                     0.,
                     180.0),
         neointimaForce(0.095,
                     0.,
                     0.,
                     0.,
                     0.,
                     68.),
         threeLayerForce(&zeroBinForce),
         solverName(solverName_)
    {
        setupThreeLayerForce();
        logger() << "Solver \"" << solverName << "\" created at " << Timer::printTime() << std::endl;

    }

    ~AgentSolver() {
    }

    void loadConfig(std::string const& configFileName) {
        configDir = parseFileName(configFileName).directory;
        dataInputDir = configDir;
        outputDir = createFileName(configDir, std::string("data.") + solverName, "");
        config().readFile(configFileName);
    }

    void setDataInputDir (std::string const& dataInputDir_) {
        dataInputDir = dataInputDir_;
    }

    void setInFileName(std::string const& inFileName_) {
        inFileName = inFileName_;
    }

    void setOutFileName(std::string const& outFileName_) {
        outFileName = outFileName_;
    }

    void setBaseNameFromIn() {
        baseName = parseFileName(inFileName).baseName;
    }

    void setBaseNameFromOut() {
        baseName = parseFileName(outFileName).baseName;
    }

    void readIOFileNames(std::string const& inFileNameParam_, std::string const& outFileNameParam_="") {
        // agent input file
        setInFileName(config().getValue<std::string>(inFileNameParam_));
        if (outFileNameParam_ != "") {
            setOutFileName(config().getValue<std::string>(outFileNameParam_));
            setBaseNameFromOut();
        } else {
            setBaseNameFromIn();
        }
    }

    std::string getBaseName() {
        return baseName;
    }

    std::string getInFileName() {
        return inFileName;
    }

    std::string getOutFileName() {
        return outFileName;
    }


    void readCharLengthScales() {
        // characteristic length scale
        setCharLengthScales(config().getValue<double>("smc_mean_rad"),
                            config().getValue<double>("smc_rad_sigma"));
    }

    void setCharLengthScales(double L_, double sigma_) {
        L = L_;
        sigma = sigma_;
    }

    // load agents from file and put them in AgentContainer
    void readCsv3DFile(std::string const& filename) {
        AgentFileReader<3>::readCsvFile(createFileName(dataInputDir, filename), AgentFactory3D(), agents);
        nAgents = agents.count();
        logger() << nAgents << " agents loaded from csv." << std::endl;
        if(nAgents == 0)
            logger() << "WARNING: NO AGENTS LOADED" << std::endl;
    }

    void readDat3DFile(std::string const& filename) {
        logger() << "Reading initial cells from file " << dataInputDir << "/" << filename << std::endl;

        AgentFileReader<3>::readFile(createFileName(dataInputDir, filename), AgentFactory3D(), agents);
        nAgents = agents.count();
        logger() << nAgents << " agents loaded from dat." << std::endl;
        logger() << "Including: " << agents.count(tSMC3D) << " SMCs, "
                                     << agents.count(tIEL3D) << " IEL agents, "
                                     << agents.count(tEEL3D) << " EEL agents." << std::endl;
        if(nAgents == 0)
            logger() << "WARNING: NO AGENTS LOADED" << std::endl;
    }

    /// load stent from file
    void loadGeneratedStent() {
        const std::string stentFileName = config().getValue<std::string>("deploy_stent_file");
        stent = std::make_unique< SimpleStent3D >(createFileName(dataInputDir, stentFileName), 0.0, 0.0 );
        stent->moveToCenterX((agents.getMaxCoord(0) + agents.getMinCoord(0)) / 2.);
        stent->exportObstacles(agents);
        logger() << stent->count() << " stent obstacles loaded" << std::endl;
    }

    /// load balloon from file
    void loadGeneratedBalloon() {
        const std::string balloonFileName = config().getValue<std::string>("deploy_balloon_file");
        balloon = std::make_unique< SimpleBalloon3D >(createFileName(dataInputDir, balloonFileName), 0.0, 0.0 );
        balloon->moveToCenterX((agents.getMaxCoord(0) + agents.getMinCoord(0)) / 2.);
        balloon->exportObstacles(agents);
        logger() << balloon->count() << " balloon obstacles loaded" << std::endl;
    }

    void removeBalloon() {
        logger() << "Removing balloon" << std::endl;
        balloon->remove(agents);
    }

    // range should be enough to cover all configurations
    void initStaticNeighbourDetector() {
        nbDetector = std::make_unique< StaticNeighbourDetector<agent_t> >(2.7 * L + 2 * sigma);
        nbDetector->init(agents.getAgentVector());
    }

    void initVerletNeighbourDetector() {
        nbDetector = std::make_unique< VerletNeighbourDetector<agent_t> >(2 * L, 1.5*2*L); //works in general case, large neighbourhoods
    }

    void makeBondsFromNeighbours() {
        /// a bit dirty, can be a container-wide thing or pre-set in csv
        nAgents = agents.count();
        for (size_t iAgent = 0; iAgent < nAgents; iAgent++) {
            //agents.getAgentVector()[iAgent]->addNeighboursToBonds();
            agents.getAgentVector()[iAgent]->addNeighboursToBonds(true, 12); // max bonds set to 12 since that's the number of neighbours in a closest packing
        }
        logger() << "Bonds initialized from neighbours" << std::endl;
    }

    /// Detect surface agents. Required for any other operations with the surface
    void detectSurface() {
        // SURFACE DETECTION
        Timer timer;
        logger() << std::endl << "Starting surface detection..." << std::endl;
        timer.start();
        NeighbourCache<agent_t>* nbCacheForCiRule = new NeighbourCache<agent_t > (agents.getAgentVector());
        FDetectSurfaceFromNeighbours<agent_t> surfaceDetectorFunctor(L * 7, 150, &nbCacheForCiRule); //200 for the denser (and deprecated) hexahedral lattice
        agents.forAll(surfaceDetectorFunctor);
        timer.stop();
        logger() << "Surface detection complete, duration: " << timer.printElapsedTime() << std::endl;
        delete nbCacheForCiRule;
    }

    void countSurface () {
        size_t surfaceCount = 0;
        /// Count surface, set IEL and EEL
        for (size_t iAgent = 0; iAgent < agents.count(); iAgent++) {
            agent_t* agent = agents.getAgentVector()[iAgent];
            if(agent->isSurface()){
                surfaceCount++;
            }
        }
        logger() << surfaceCount << " surface agents" << std::endl;
    }

    void surfaceToMembranes() {
        size_t surfaceCount = 0;
        /// Count surface, set IEL and EEL
        for (size_t iAgent = 0; iAgent < agents.count(); iAgent++) {
            agent_t* agent = agents.getAgentVector()[iAgent];
            if(agent->isSurface()){
                surfaceCount++;

                if(centerline.isFacingCenterline(agent->getPos(), *(agent->getSurfaceNormal()))) {
                    //agent->setAffected(false);
                    agents.replace(iAgent, new IEL3D(agent->getPos(), agent->getPos0(), agent->getR()));
                }

                else {
                    agent->unsetSurfaceNormal();
                    agents.replace(iAgent, new EEL3D(agent->getPos(), agent->getPos0(), agent->getR()));
                }
            }
        }
        logger() << surfaceCount << " surface agents" << endl;
    }

    /// Set inner surface immobile
    void setSurfaceImmobile() {
        logger() << "Fixing the inner surface in place" << std::endl;
        for (size_t iAgent = 0; iAgent < agents.count(); iAgent++) {
            agent_t* agent = agents.getAgentVector()[iAgent];
            if(agent->isSurface()){
                if(centerline.isFacingCenterline(agent->getPos(), *(agent->getSurfaceNormal()))) {
                    agent->setMobility(point_t(0.,0.,0.));
                }
            }
        }
    }

    void readBondsFromFile(std::string const& filename) {
        const std::string bondFileName = createFileName(dataInputDir, filename);
        NBReader<nDim>::readFile(bondFileName, agents, true);
    }

    void readCenterlineFromFile() {
        const std::string centerlineFileName = config().getValue<std::string>("centerline_file");
        readCenterlineFromFile(centerlineFileName);
    }

    void readCenterlineFromFile(std::string const& filename) {
        centerline.readFromFile(createFileName(dataInputDir,filename));
    }

    void readWallDisplacements() {
        logger() << "Loading vessel wall displacements..." << std::endl;
        const std::string displacementsFileName =config().getValue<std::string>("displacements_file");
        displacements.readFromFile(createFileName(dataInputDir,displacementsFileName));
        displacements.printCount();
    }

    void readStentWithDisplacements(std::string const& stentFilenameInitParam, std::string const& stentFilenameFinalParam, std::string const& displacementsFilenameParam) {
        const std::string stentFileName              = config().getValue<std::string>(stentFilenameInitParam);
        const std::string finalStentFileName         = config().getValue<std::string>(stentFilenameFinalParam);
        const std::string stentDisplacementsFileName = config().getValue<std::string>(displacementsFilenameParam);
        logger() << "Loading stent displacements..." << std::endl;
        displacementsStent.readFromFile(createFileName(dataInputDir,stentDisplacementsFileName));
        displacementsStent.printCount();

        // load stent (initial configuration)
        logger() << "Loading initial stent configuration..." << std::endl;
        stent = std::make_unique< SimpleStent3D >(createFileName(dataInputDir, stentFileName), 0.0, 0.0 );
        stent->exportObstacles(agents);
        logger() << stent->count() << " obstacles loaded" << std::endl;

        // also read the final configuration for later
        finalStent = std::make_unique< SimpleStent3D >(createFileName(dataInputDir, finalStentFileName), 0.0, 0.0 );
    }

    // set initial positions of cells to current positions
    void resetPos0() {
        agents.forAll(FResetPos0<agent_t>() );
    }

    void trimBoundaries(double trimDistance) {
        agents.remove(PSelectCloseToBoundary<agent_t>(0, agents.getMinCoord(0), agents.getMaxCoord(0), trimDistance));
        logger() << "After trimming longitudinal boundaries: " << agents.count() << " agents left." << std::endl;
        logger() << "Including: " << agents.count(tSMC3D) << " SMCs, "
                                     << agents.count(tIEL3D) << " IEL agents, "
                                     << agents.count(tEEL3D) << " EEL agents." << std::endl;
   }

    void setBoundaryMobility(std::string const& mobility_x, std::string const& mobility_y, std::string const& mobility_z, std::string const& immobility_range_par = "") {
        const double mX = config().getValue<double>(mobility_x);
        const double mY = config().getValue<double>(mobility_y);
        const double mZ = config().getValue<double>(mobility_z);
        const double immobility_range = config().hasValue(immobility_range_par) ? config().getValue<double>(immobility_range_par) : L*4;
        setBoundaryMobility(mX, mY, mZ, immobility_range);
    }

    void setBoundaryMobility(const double mX, const double mY, const double mZ, const double immobility_range) {
            agents.forAll(PSelectCloseToBoundary<agent_t>(0, agents.getMinCoord(0), agents.getMaxCoord(0), immobility_range),
                          FSetMobility<agent_t>(point_t(mX, mY, mZ) ) );
    }

    void readIterationParameters(std::string const& maxIter_, std::string const& vtkIter_ = "", std::string const& datIter_ = "") {
        maxIter     = config().getValue<int>(maxIter_);
        vtkIter     = config().hasValue(vtkIter_) ? config().getValue<int>(vtkIter_) : 0;
        datIter     = config().hasValue(datIter_) ? config().getValue<int>(datIter_) : 0;
    }

    void setIterationParametersFromDisplacements(std::string const& vtkIter_ = "", std::string const& datIter_ = "") {
        if (displacements.stepsCount() != displacementsStent.stepsCount()) {
            if(displacementsStent.stepsCount() == 0) {
                logger() << "No stent displacements, basing the number of steps on the artery displacements" << std::endl;
            }
            else {
                logger() << "WARNING: Step counts for the artery and for the stent are different (" << displacements.stepsCount() << " and " << displacementsStent.stepsCount() << "), using the larger value" << std::endl;
            }
        }
        /// iteration parameters
        maxIter     = std::max(displacements.stepsCount(), displacementsStent.stepsCount()) + 2;
        vtkIter     = config().hasValue(vtkIter_) ? config().getValue<int>(vtkIter_) : 0;
        datIter     = config().hasValue(datIter_) ? config().getValue<int>(datIter_) : 0;
    }

    void setVtkIter(const double vtkIter_) {
        vtkIter = vtkIter_;
    }

    void setDatIter(const double datIter_) {
        datIter = datIter_;
    }

    /// Output values for .vtp files; read from absmc config
    void readVtpParameters(std::string const& scalarsParam, std::string const& vectorsParam) {
        setVtpParameters(config().getValue<std::string>(scalarsParam), config().getValue<std::string>(vectorsParam));
    }

    /// Output values for .vtp files; directly set to input strings
    void setVtpParameters(std::string const& scalars, std::string const& vectors) {
        vtpScalars = scalars;
        vtpVectors = vectors;
    }

    void writeAllOutputs(std::string const& prefix) {
        writeDatOutput(prefix);
        writeVtkOutput(prefix);
        logger() << "Output from " << solverName << " successfully written" << std::endl;
    }

    void writeDatOutput(std::string const& prefix) {
        agents.writeFile(createFileName(configDir, prefix+baseName, "dat") );
        NBWriter<nDim>::writeFile(agents, createFileName(configDir, prefix+baseName+"_nb", "dat"));
    }

    void writeDatOutputNoBonds(std::string const& prefix) {
        agents.writeFile(createFileName(configDir, prefix+baseName, "dat") );
    }

    void writeVtkOutput(std::string const& prefix) {
        graphics::VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(), createFileName(configDir, prefix+baseName, "vtp"), vtpScalars, vtpVectors);
    }

    /// Shuffle agent positions; the direction of the displacement vector is uniformly distributed,
    /// and its norm is Gauss-distributed with specified mean and sigma. Only the specified fraction of agents is shuffled,
    /// the remainder is left at their initial position and set immobile. This allows to fix a certain number of dof when
    /// studying if the ensemble returns to a pre-shuffling equilibrium configuration.
    /// Agents already marked as immobile are not displaced either; note that if their number is non-zero,
    /// the fraction of agents actually displaced will differ from the specified fraction.
    void shuffleAgents(double meanShift, double sigmaShift, double shuffleFraction)
    {
        auto agentVector = agents.getAgentVector();
        nAgents = agentVector.size();
        size_t nFixed = 0;
        for (size_t iAgent=0; iAgent<nAgents; iAgent++) {

            const double shift = util::Random::getNextGauss(meanShift, sigmaShift);
            const double theta = util::Random::getNext(math::pi);
            const double phi   = util::Random::getNext(2.0*math::pi);

            point_t newPos = agentVector[iAgent]->getPos();

            newPos[0] += shift * sin(theta) * cos(phi);
            newPos[1] += shift * sin(theta) * sin(phi);
            newPos[2] += shift * cos(theta);

            // displace only the specified fraction of agents, fix the remainder.
            // do not displace agents which were already fixed.
            const double p = util::Random::getNext(1.0);
            if (p < shuffleFraction && agentVector[iAgent]->getMobility() == point_t(1.0, 1.0, 1.0) ) {
                agentVector[iAgent]->setPos(newPos);
            }
            else {
                agentVector[iAgent]->setMobility(point_t(0.0, 0.0, 0.0) );
                nFixed++;
            }
        }
        logger() << (nAgents-nFixed) << " agents displaced randomly with meanShift=" << meanShift << ", sigmaShift=" << sigmaShift << ". "
                  << nFixed << " agents undisplaced and set immobile." << std::endl;
    }

    void shuffleAgentsByConfig() {
        const double meanShift  = config().getValue<double>("eq_mean_shift");
        const double sigmaShift = config().getValue<double>("eq_sigma_shift");
        const double shuffleFraction = config().getValue<double>("eq_shuffle_fraction");
        shuffleAgents(meanShift, sigmaShift, shuffleFraction);
    }

    Poly6thPowerAttractionForce3D createPolynomialForce(std::string const& c1_par,
                                                        std::string const& c2_par,
                                                        std::string const& c3_par,
                                                        std::string const& c4_par,
                                                        std::string const& c5_par,
                                                        std::string const& c6_par) {
        const double c1 = config().getValue<double>(c1_par);
        const double c2 = config().getValue<double>(c2_par);
        const double c3 = config().getValue<double>(c3_par);
        const double c4 = config().getValue<double>(c4_par);
        const double c5 = config().getValue<double>(c5_par);
        const double c6 = config().getValue<double>(c6_par);
        return Poly6thPowerAttractionForce3D(c1, c2, c3, c4, c5, c6);
    }

    UnaryTrajectoryForce createWallTrajectoryForce() {
        const double iterTime = config().getValue<double>("deploy_iter_time");
        double unaryMagnitudeTrajForce = (1. / iterTime);
        unaryMagnitude = std::max(unaryMagnitude, unaryMagnitudeTrajForce / 10);
        return UnaryTrajectoryForce (unaryMagnitudeTrajForce, displacements);
    }

    UnaryTrajectoryForce createStentTrajectoryForce() {
        const double iterTime = config().getValue<double>("deploy_iter_time");
        double unaryMagnitudeTrajForce = (1. / iterTime);
        unaryMagnitude = std::max(unaryMagnitude, unaryMagnitudeTrajForce / 10);
        // only take the closest trajectory point into account
        return UnaryTrajectoryForce (unaryMagnitudeTrajForce, displacementsStent, 1);
    }

    UnaryRadialDisplacementForce createRadialStentDisplacementForce() {
        const double iterTime = config().getValue<double>("deploy_iter_time");
        const double deploymentDepth = config().getValue<double>("deployment_depth");
        const double dr = (deploymentDepth) / maxIter;
        logger() << "Stent dr per macro iter = " << dr << std::endl;

        unaryMagnitude = dr / iterTime; // stent displacement
        return UnaryRadialDisplacementForce (unaryMagnitude, std::make_shared<Centerline<3> >(centerline));
    }

    void setForces(BinaryForce< AgentBase<nDim> > * binForce_,
                   BondForce< AgentBase<nDim> > * bondForce_,
                   UnaryForce < AgentBase<nDim> > * uForce_) {
        binForce = binForce_;
        bondForce = bondForce_;
        uForce = uForce_;
    }

    void calculateCharForce() {
        /// system-characteristic force
        const double compressionMagnitude = std::abs(binForce->calculateForce(L, L, L));
        const double stretchMagnitude = std::abs(bondForce->calculateForce(L, L, 2.3 * L));
        binaryMagnitude = std::max(compressionMagnitude, stretchMagnitude);
        charForce = std::max(binaryMagnitude, unaryMagnitude); //REAL UNITS

        // maximal displacement
        solverMaxDispl = 0.05*L; //REAL UNITS //0.03L?
        // maximal time step (from rough stability estimation)
        const double maxDt = solverMaxDispl / charForce;
        // maximal time step passed to solver
        solverMaxDt = 500. * maxDt;

        logger() << "charForce         = " << charForce << std::endl;
        logger() << "maxDt             = " << maxDt << std::endl;
        logger() << "solverMaxDt       = " << solverMaxDt << std::endl;
        logger() << "solverMaxDispl    = " << solverMaxDispl << std::endl;

        logger() << std::endl << "Interaction and bond forces for various distances" << std::endl;
        for (double dist = 1.0; dist < 3.0; dist += 0.1) {
            logger() << "Distance = " << dist << " R; Interaction force = " << binForce->calculateForce(L, L, L*dist) <<
                                                    "; Bond force = " << bondForce->calculateForce(L, L, L*dist) << std::endl;
        }
    }

    void initializeIntegrator() {
        std::ofstream integratorLogFile(createFileName(outputDir, "fr", "dat").c_str() );
        integrator = std::make_unique< RungeKuttaIntegratorMPI<agent_t> >(*uForce, *binForce, *bondForce, L, solverMaxDispl, solverMaxDt, *nbDetector, integratorLogFile);
    }

    void setForceResidualMaxNormController(std::string const& epsParam) {
        // desired convergence level
        const double eps      = config().getValue<double>(epsParam);
        setForceResidualMaxNormController(eps);
    }

    void setForceResidualMaxNormController(const double eps) {
        logger() << "eps               = " << eps << std::endl;
        // can set 1st param to limit number of iterations, e.g. 1000
        controller = std::make_unique< ForceResidualMaxNormController >(std::numeric_limits<double>::infinity(), charForce, eps);
    }

    void setFixedIntervalController(std::string const& iterTimeParam) {
        const double iterTime = config().getValue<double>(iterTimeParam);
        controller = std::make_unique< FixedIntervalController >(iterTime);
    }

    void setDeploymentRules() {
        /// set rules for SMC agents (necrosis rule)
        const double smcMaxStress = config().getValue<double>("smc_max_stress");
        smcNecrosisRule = std::make_shared< SMCNecrosisRule3D >(smcMaxStress);

        smcRule = std::make_shared< CompositeRule<SMC3D> > ();
        smcRule->add(smcNecrosisRule);

        agents.forAll(tSMC3D, FSetAgentRule<SMC3D>(smcRule) );

        /// set rules for IEL agents (iel breaking rule)
        const double ielMaxStrain = config().getValue<double>("iel_max_strain");
        const double ielMaxStress = config().getValue<double>("iel_max_stress");
        ielBreakingRule = std::make_shared< IELBreakingRule3D >(ielMaxStrain , ielMaxStress);
        agents.forAll(tIEL3D, FSetAgentRule<IEL3D>(ielBreakingRule) );
    }

    /// strain along X axis
    void applyEngineeringStrain(double cur_strain) {
        agents.forAll(FApplyEngineeringStrain<agent_t>(cur_strain, 0));
    }

    /// Stress for a cubic sample spretched along X
    double calculateStress(const double axis = 0) {
        const double minX = agents.getMinCoord(0);
        const double minY = agents.getMinCoord(1);
        const double minZ = agents.getMinCoord(2);
        const double maxX = agents.getMaxCoord(0);
        const double maxY = agents.getMaxCoord(1);
        const double maxZ = agents.getMaxCoord(2);

        const double crossSectionArea = (maxY - minY) * (maxZ - minZ);

        auto agentVector = agents.getAgentVector();
        nAgents = agentVector.size();
        double forceSumLow = 0.;
        double forceSumHigh = 0.;
        for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
            point_t pos = agentVector[iAgent]->getPos();
            if ((pos[0] < (minX + maxX) / 2) && (!agentVector[iAgent]->getMobile())) {
                forceSumLow += agentVector[iAgent]->getForce()[0];
            }
            if ((pos[0] > (minX + maxX) / 2) && (!agentVector[iAgent]->getMobile())) {
                forceSumHigh -= agentVector[iAgent]->getForce()[0];
            }
        }
        double stress = (forceSumLow + forceSumHigh) / 2. / crossSectionArea;
        logger() << "Total force low boundary " << forceSumLow << " N, high boundary " << forceSumHigh << std::endl;
        logger() << "Cross-sectional stress is " << stress << " MPa" << std::endl;
        return stress;
    }

    void singleIteration(size_t iter, Timer& tIter, Timer& tTotal) {
        if (vtkIter!=0 && iter%vtkIter==0) {
            tIter.stop();
            tTotal.stop();
            graphics::VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(),
                                                           createFileName(outputDir, baseName, "vtp", iter, 6),
                                                           vtpScalars, vtpVectors);
            logger() << "iter=" << iter << ",  t(" << vtkIter << " it)="
                      << tIter.printElapsedTime() << ",  t(total)=" << tTotal.printElapsedTime()
                      << ", nSMC=" << agents.count(tSMC3D) << ", nIEL=" << agents.count(tIEL3D) << std::endl;
            tIter.start();
        }
        if (datIter!=0 && iter%datIter==0) {
            agents.writeFile(createFileName(outputDir, baseName, "dat", iter, 6) );
            NBWriter<nDim>::writeFile(agents, createFileName(outputDir, baseName+"_nb", "dat", iter, 6));
        }
        integrator->integrate(agents.getAgentVector(), *controller);
        /// exec agent rules (iel breaking, smc necrosis)
        agents.execAgentRules();
    }

    // iterates until maxIter is reached
    void iterate() {
        Timer tIter;
        Timer tTotal;
        for (size_t iter=0; iter<maxIter; iter++) {
            singleIteration(iter, tIter, tTotal);
        }
        tTotal.stop();
        logger() << "Iterations complete, total elapsed time " << tTotal.printElapsedTime() << std::endl;
    }

    void updateTrajectories(size_t step, UnaryTrajectoryForce& trajForce, UnaryTrajectoryForce& trajForceStent) {
        trajForce.updateAllForces(agents.getAgentVector(), step);
        trajForceStent.updateAllForces(agents.getAgentVector(), step);
    }

    // integrate and keep the time-dependent forces up to date.
    void iterateWithTrajectories(UnaryTrajectoryForce& trajForce) {
        Timer tIter;
        Timer tTotal;
        trajForce.updateAllForces(agents.getAgentVector(), 0);

        for (size_t iter=0; iter<maxIter; iter++) {
            singleIteration(iter, tIter, tTotal);
            trajForce.updateAllForces(agents.getAgentVector(), iter + 1);
        }
    }

    /// set trajectory forces to.zero and do a final iteration
    void finalIterateWithTrajectories(UnaryTrajectoryForce& trajForce) {
        Timer tIter;
        Timer tTotal;
        trajForce.updateAllForces(agents.getAgentVector(), maxIter + 1);
        singleIteration(maxIter, tIter, tTotal);
    }


    // integrate and keep the time-dependent forces up to date, also for the stent.
    void iterateWithTrajectories(UnaryTrajectoryForce& trajForce, UnaryTrajectoryForce& trajForceStent) {
        Timer tIter;
        Timer tTotal;
        updateTrajectories(0, trajForce, trajForceStent);

        for (size_t iter=0; iter<maxIter; iter++) {
            singleIteration(iter, tIter, tTotal);

            /// add the stent when the trajectory has ended
            /// switch the deformed stent we got from the initial state with the final state from FE
            if(iter + 1 == displacements.stepsCount()) {
                stent->remove(agents);
                // load stent (final configuration)
                logger() << "Loading final stent configuration..." << std::endl;
                finalStent->exportObstacles(agents);
                logger() << finalStent->count() << " obstacles loaded" << std::endl;
            }
            updateTrajectories(iter + 1, trajForce, trajForceStent);
        }
    }

    /// This function adds fenestrations to the IEL, replacing IEL agents with SMC agents with a set [probability]
    void makeFenestrations (double probability) {
        std::vector<agent_t*> & agentsVector = agents.getAgentVector();
        nAgents = agentsVector.size();
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0., 1.);
        for (size_t iAgent = 0; iAgent < nAgents; iAgent++) {
            if (agentsVector[iAgent]->getTypeId() == tIEL3D) {
                /// Change to SMC with a set probability
                if (distribution(generator) <= probability) {
                    agents.replace(iAgent, new SMC3D(agentsVector[iAgent]->getPos(), agentsVector[iAgent]->getPos0(), agentsVector[iAgent]->getR()));
                }
            }
        }
    }

    void setEelImmobile() {
        agents.forAll(tEEL3D, FSetMobility<agent_t>(point_t(0.0, 0.0, 0.0) ) );
    }

    void setBalloonDenudation(double balloonProtrusion) {
        //set endothelium to denuded if the cell is under the balloon; set the rest inactive
        logger() << "Denuding endothelium under the balloon..." << std::endl;
        const double stentMin = agents.getMinCoord(0, tObstacle3D);
        const double stentMax = agents.getMaxCoord(0, tObstacle3D);
        logger() << "Stent from " << stentMin << " to " << stentMax << std::endl;
        agents.forAll(PSelectByDisplacement<agent_t > (2.*L), FDenudeSMC3DEndothelium());
        //agents.forAll(PSelectInterval<agent_t > (0, stentMin - balloonProtrusion, stentMax + balloonProtrusion), FDenudeSMC3DEndothelium());
    }

    /// set SMC agents outside of balloon-denuded area biologically inactive
    void setSMCsInactive(double balloonProtrusion) {
        const double stentMin = agents.getMinCoord(0, tObstacle3D);
        const double stentMax = agents.getMaxCoord(0, tObstacle3D);
        agents.forAll(PSelectOutsideInterval<agent_t > (0, stentMin - balloonProtrusion, stentMax + balloonProtrusion), FSetSMC3DInactive());
        //agents.forAll(PSelectCloseToBoundary<agent_t > (0, agents.getMinCoord(0), agents.getMaxCoord(0), BR), FSetSMC3DInactive());
    }

    /// Calculate how much of the surface should be EC-covered at this timestep
    double calculateEndothelialPercentage(const double timeSinceStart, const double endPoint) {
        std::vector <double> timePoints {0., 24. * 3., 24. * endPoint, std::numeric_limits<double>::max() / 2.};
        std::vector <double> coverage   {0., 0.59    , 1.0,            1.0};

        size_t curIndex = 0;
        for(; timePoints[curIndex] <= timeSinceStart; curIndex++);

        const double endothelialCoverage = coverage[curIndex - 1] +
                                           (coverage[curIndex] - coverage[curIndex - 1]) *
                                           (timeSinceStart - timePoints[curIndex - 1]) / (timePoints[curIndex] - timePoints[curIndex - 1]);

        return endothelialCoverage;
    }


protected:
    AgentContainer<nDim> agents;
    size_t nAgents;

    Centerline<nDim> centerline;
    TrajectorySet<nDim> displacements;
    TrajectorySet<nDim> displacementsStent;


    std::unique_ptr< SimpleStent3D > stent;
    std::unique_ptr< SimpleStent3D > finalStent;
    std::unique_ptr< SimpleBalloon3D > balloon;
    std::unique_ptr< NeighbourDetector<agent_t> > nbDetector;

    const std::string solverName;
    std::string configDir, dataInputDir, outputDir;
    std::string inFileName, outFileName, baseName, vtpScalars, vtpVectors;
    double L, sigma;

    BinaryForce< AgentBase<nDim> > * binForce;
    BondForce< AgentBase<nDim> > * bondForce;
    UnaryForce < AgentBase<nDim> > * uForce;
    std::unique_ptr< RungeKuttaIntegratorMPI<agent_t> > integrator;
    std::unique_ptr< IntegratorController > controller;

    size_t maxIter, vtkIter, datIter;

    double charForce, binaryMagnitude, unaryMagnitude, solverMaxDispl, solverMaxDt;

    std::shared_ptr< SMCNecrosisRule3D > smcNecrosisRule;
    std::shared_ptr< SMCContactInhibitionRule3D > smcContactInhibitionRule;
    std::shared_ptr< SMCNitricOxideRule3D > smcNitricOxideRule3D;
    std::shared_ptr< SMCCellCycleRule3D > smcCellCycleRule;
    std::shared_ptr< ECMProductionRule3D > ecmProductionRule3D;
    std::shared_ptr< SMCBondBreakingRule3D > smcBondBreakingRule;

    std::shared_ptr< CompositeRule<SMC3D> > smcRule;
    std::shared_ptr< IELBreakingRule3D > ielBreakingRule;


};

} // end namespace absmc

#endif
