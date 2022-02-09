#ifndef FLOW_SOLVER_H
#define FLOW_SOLVER_H

#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <cstdlib>
#include "palabos3D.h"
#include "palabos3D.hh"
#include <libmuscle/libmuscle.hpp>
#include <ymmsl/ymmsl.hpp>
#include <util/logger.h>

using namespace std;
typedef double T;
using util::logger;
using namespace plb;
using namespace plb::util;

#define DESCRIPTOR descriptors::D3Q19Descriptor

enum VesselTypeID { tAny, tRectangular, tCylindrical, tCurved };

#include "initialization.h"
#include "flowFileInput.h"
#include "flowFileOutput.h"


using libmuscle::Data;
using libmuscle::DataConstRef;
using libmuscle::Instance;
using libmuscle::Message;
using ymmsl::Operator;

class FlowSolver {
public:
    FlowSolver(Instance& instance_);
    ~FlowSolver();

    void receiveBoxSize();
    void receiveMask();

    void initialize();
    void boundaryConditionSetup();
    void generateParameters(const T uConst);
    void execute(int argc, char* argv[]);
    void singleIteration();

    void writeOutput(bool writeToFile);
    void sendShearStress(MultiScalarField3D<T>& realwsstress);

    T shearStressConversionFactor();
    T velocityConversionFactor();
    T pressureConversionFactor();
    T calculateReFromReal();

private:
    Instance* instance;
    plint N, Nx, Ny, Nz;

    double convergence_limit;
    double velocity_real;
    double dx_real;
    int64_t vtkIter;
    int64_t datFieldsIter;
    bool writeGifs;

    const double blood_nu = 2.7e-6; //blood viscosity at 36 deg. C  in [m^2 / s]
    double timeInHours;
    size_t sz, initFluidNodes, macroIter; //number of lattice cells and fluid nodes

    IsoThermalBulkDynamics<T, DESCRIPTOR>* defaultDynamics;
    MultiBlockLattice3D<T, DESCRIPTOR>* lattice;
    IncomprFlowParam<T>* parameters;
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition;

    std::vector<int32_t> geometry;          //vessel geometry as received from distributor
    MultiScalarField3D<int32_t>* charMask;  //3D mask made from the received geometry
    bool maskReceived;      // whether a new mask was received on this iteration

    MultiTensorField3D<T, 6 >* wsstress;
    VesselTypeID vesselTypeID;
};


FlowSolver::FlowSolver(Instance& instance_) :
instance(&instance_), sz(0), initFluidNodes(0), macroIter(0) {
    convergence_limit = instance->get_setting_as<double>("convergence_limit");
    velocity_real = instance->get_setting_as<double>("flow_velocity");
    dx_real = instance->get_setting_as<double>("flowdx");
    vtkIter = instance->get_setting_as<int64_t>("bf_vtk_iter");
    datFieldsIter = instance->get_setting_as<int64_t>("bf_dat_iter");
    writeGifs = instance->get_setting_as<int64_t>("write_gifs");
}


FlowSolver::~FlowSolver() {
    // Manually free objects we own.
    delete lattice;
    delete charMask;
    delete wsstress;
    delete parameters;
    delete boundaryCondition;
}

void FlowSolver::receiveMask() {
    auto geometryIn = instance->receive("geometryIn");//, dummyWSSMessage);
    size_t newSz;
    size_t nDim = 3;
    std::vector<size_t> boxSize(nDim);

    if (global::mpi().isMainProcessor()) {
        timeInHours = geometryIn.timestamp();
        logger() << "Geometry received for timepoint " << timeInHours << " hours" << std::endl;
        boxSize = geometryIn.data().shape();

        newSz = geometryIn.data().size();
        geometry.resize(newSz);
        auto obsPtr = geometryIn.data().elements<int32_t>();
        geometry.assign(obsPtr, obsPtr + newSz);

    }

    global::mpi().bCast(&timeInHours, 1);
    global::mpi().bCast(&newSz, 1);
    if (!global::mpi().isMainProcessor()) {
        geometry.resize(newSz);
    }
    global::mpi().bCast(boxSize.data(), nDim);
    global::mpi().bCast(geometry.data(), newSz);

    pcout << "Broadcasts complete" << std::endl;

    Nx = boxSize[0];
    Ny = boxSize[1];
    Nz = boxSize[2];

    /// TODO: add check for changing box size during runtime

    size_t fluidNew = 0;
    if (global::mpi().isMainProcessor())
    {
        pcout << "Counting fluid nodes" << std::endl;

        if (sz == 0) {
            sz = newSz;
        }
        else if (newSz > 0 && newSz != sz) {
            logger() << "Previous number of lattice cells: " << sz << ", new number of cells: " << newSz << std::endl;
            throw runtime_error("Can not change the number of lattice cells during simulation");
        }

        // count fluid cells, used in generateParameters for calibration
        for (size_t i = 0; i < newSz; i++) {
            if (geometry[i] == 0) {
                fluidNew++;
            }
        }
        logger() << "Received " << fluidNew << " fluid cells" << std::endl;
    }
    global::mpi().bCast(&fluidNew, 1);
    if (fluidNew == 0) {
        maskReceived = false;
    }

    // Initialize the number of fluid nodes and characteristic length
    if (initFluidNodes == 0) {
        initFluidNodes = fluidNew;

        //characteristic length (vessel width if the vessel was square)
        //Initial value is used for the scale conversion throughout the simulation
        N = (plint) sqrt(initFluidNodes / Nx);
    }

    maskReceived = true;
    pcout << "Mask receive complete" << std::endl;
}


void FlowSolver::initialize() {
    // omega is 1/tau, where tau is the relaxation time of the LB fluid
    // omega = 1/(3*nu+0.5) = 1/(3*uMax*N/Re + 0.5) = 1/(3*(2/N)*N/Re + 1/2) = 1/(6/120+1/2) = 1.82
    defaultDynamics = new IncBGKdynamics<T, DESCRIPTOR>(1/(6/calculateReFromReal() + 1/2));
    // Lattice that the computation will run on
    lattice = new MultiBlockLattice3D<T, DESCRIPTOR> (Nx, Ny, Nz, defaultDynamics);
    pcout << "Lattice is created" << std::endl;

    boundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    pcout << "Boundary conditions are initialized" << std::endl;

    lattice->toggleInternalStatistics(true);
    charMask = new MultiScalarField3D<int32_t>(Nx, Ny, Nz);
    wsstress = new MultiTensorField3D<T, 6 >(Nx, Ny, Nz);
    pcout << "Fields are initialized" << std::endl;

    std::string vesselTypeName = instance->get_setting_as<std::string>("vessel_type");
    vesselTypeID = tAny;
    if (vesselTypeName == "rectangular") vesselTypeID = tRectangular;
    if (vesselTypeName == "cylindrical") vesselTypeID = tCylindrical;
    if (vesselTypeName == "curved") vesselTypeID = tCurved;
    if (vesselTypeID == tAny) {
        throw runtime_error("Unknown vessel type. Aborting.");
    }
    pcout << "Solver initialization complete" << std::endl;
}


// Calculates Re from flow parameters in physical units
// Re is based on initial velocity, not on current one
// Calculation is done in SI units (m, s)
T FlowSolver::calculateReFromReal() {
    //Re = u_real * l_real / nu_real

    //corresponds to unit length in dimensionless system (N cells), characteristic size of the system
    T x_characteristic = dx_real * 0.001 * N;
    T Re = velocity_real * x_characteristic / blood_nu;
    pcout << "Reynolds number calculated as " << Re << std::endl;
    return Re;
}


//Calculates WSS in real units (Pa) from strain rate in lattice units
T FlowSolver::shearStressConversionFactor() {
    //wss = nu * rho * du / dy
    T rho_real = 1000.; //blood is approximated with water
    T u_conversion_factor = velocity_real / parameters->getLatticeU();

    //corresponds to unit length (one cell) in LB system
    T x_conversion_factor = dx_real * 0.001;
    T nu_real = blood_nu;

    T conversion_factor = nu_real * rho_real  * u_conversion_factor / x_conversion_factor;
    return conversion_factor;
}


//Converts velocity from lattice units to real units (m/s)
T FlowSolver::velocityConversionFactor() {
    T u_conversion_factor = velocity_real / parameters->getLatticeU();
    return u_conversion_factor;
}


/// Converts pressure from lattice units (rho, density) to real units (Pa)
/// Incompressible formulation only
T FlowSolver::pressureConversionFactor() {

    T rho_to_P_conversion_factor = DESCRIPTOR<T>::cs2;
    T rho_conversion_factor = 1000; //average lattice rho is 1; blood is approximated with water
    T u_conversion_factor = velocity_real / parameters->getLatticeU();

    /// Pa = N / m2 = kg * m / s2 / m2 = kg / m / s2 = kg / m3  *  m2 / s2
    T conversion_factor = rho_to_P_conversion_factor * rho_conversion_factor * u_conversion_factor * u_conversion_factor;
    return conversion_factor;
}


void FlowSolver::generateParameters(const T uConst)
{
    //1 (dimensionless) <=> N * dx_phys (mm)

    T Re = calculateReFromReal();
    T uMax = 0.04; // Needed to avoid compressibility errors.

    parameters = new IncomprFlowParam<T> (uMax, Re, N,
                   (T)(Nx - 1) / (T) N, // lx
                   (T)(Ny - 1) / (T) N, // ly
                   (T)(Nz - 1) / (T) N // lz
    );

    if (parameters->nStep((T) 1) == 0) {
        logger() << "The time step is 0 but should be larger, so uMax > 2*N. Given are uMax = " << uMax << ", N = " << N << std::endl; //wrong
        throw runtime_error("The time step is 0 but should be larger.");
    }
}


void FlowSolver::boundaryConditionSetup() {
    const T uMax = parameters->getLatticeU();
    const plint N = parameters->getResolution();

    Box3D inlet(0, 0, 1, Ny-2, 1, Nz-2);
    Box3D outlet(Nx-1, Nx-1, 1, Ny-2, 1, Nz-2);

     // set inlet boundary condition
    boundaryCondition->setVelocityConditionOnBlockBoundaries(*lattice, inlet, boundary::dirichlet);
    // set outlet boundary condition
    //boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, outlet, boundary::outflow);
    boundaryCondition->setPressureConditionOnBlockBoundaries(*lattice, outlet, boundary::dirichlet);
    pcout << "Inflow and outflow BCs set" << std::endl;
    //set inlet velocity profile
    setBoundaryVelocity (
        *lattice, inlet,
        PoiseuilleVelocity3D<T,DESCRIPTOR>(*parameters, *charMask) );

    const T poiseuillePressure = - 8.* uMax*uMax / parameters->getRe()
                / (N / 2) * Nx;

    const T outflowDensity = poiseuillePressure * DESCRIPTOR<T>::invCs2 + (T)1;
    pcout << "Outflow density (lattice units): " << outflowDensity << std::endl;
    setBoundaryDensity (*lattice, outlet, ConstantDensity<T>(outflowDensity));
}


void FlowSolver::sendShearStress(MultiScalarField3D<T>& realwsstress) {
    // Sending shear stress to muscle
    size_t arraysize = realwsstress.getBoundingBox().nCells() * realwsstress.sizeOfCell() / sizeof (T);
    T *mydata = new T[arraysize];
    serializerToSink(realwsstress.getBlockSerializer(realwsstress.getBoundingBox(), IndexOrdering::forward), new WriteToSerialArray<T > (mydata, (plint) arraysize));

    if (global::mpi().isMainProcessor()) {
        auto data = Data::grid(mydata, {Nx, Ny, Nz}, {"x", "y", "z"});
//        for (std::size_t i = 0; i < arraysize; ++i) {
//            data[i] = mydata[i];
//        }
        instance->send("shearStressOut", Message(timeInHours, data));
    }
    delete [] mydata;
}


void FlowSolver::writeOutput(bool writeToFile) {
    auto strainrate = computeStrainRate(*computeVelocity(*lattice));
    // computeStrainRateFromStress(*lattice, *wsstress, lattice->getBoundingBox());
    MultiScalarField3D<T> strainRateNorm = *computeSymmetricTensorNorm(*strainrate);
    //MultiScalarField3D<T> realwsstress = *extractComponent(wsstress, lattice->getBoundingBox(), 0);
    //RETURN VALUE IN Pa
    //(multiply by 10 to convert from Pa to dynes/cm^2)
    MultiScalarField3D<T> realwsstress = *multiply(/*10. * */shearStressConversionFactor(), strainRateNorm, lattice->getBoundingBox());

    pcout << "Sending shear stress at time step " << macroIter << std::endl;

    sendShearStress(realwsstress);

    if (writeToFile) {
        writeVTK(*lattice, *parameters, macroIter);
    }
    pcout << "Saving the state of the simulation" << std::endl;

    // Write data
    if ((datFieldsIter != 0) && (macroIter % datFieldsIter == 0)) {
        plb_ofstream ofile((createFileName("wssfield", macroIter , 4)+".dat").c_str());
        //plb_ofstream ofilereynolds((createFileName("reynolds", macroIter , 4)+".dat").c_str());
        //plb_ofstream ofilelspeed((createFileName("latticespeed", macroIter , 4)+".dat").c_str());
        ofile << realwsstress << endl;
        //ofilereynolds << parameters->getRe() << endl;
        //ofilelspeed << parameters->getLatticeU() << endl;

        /// write pressure (Pa)
        std::unique_ptr<MultiScalarField3D<T> > pressure_out = multiply(pressureConversionFactor(), *computeDensity(*lattice));
        plb_ofstream ofilepressure((createFileName("pressurefield", macroIter , 4)+".dat").c_str());
        ofilepressure << *pressure_out << std::endl;

        /// write velocity (m/s)
        std::unique_ptr<MultiTensorField3D<T, 3> > velocity = multiply(velocityConversionFactor(), *computeVelocity(*lattice));
        plb_ofstream ofilevelocity((createFileName("velocityfield", macroIter , 4)+".dat").c_str());
        ofilevelocity << *velocity << std::endl;
    }
}


void FlowSolver::singleIteration()
{
    // Define other parameters
    const T maxT = 1000.0; //equivalent to 1 second in real life with these parameters if I am not wrong
    const plint imSize = 600;
    global::timer("main").reset();

    pcout << "Lattice U is " << parameters->getLatticeU() << ", lattice Nu is " << parameters->getLatticeNu() << std::endl;
    pcout << "defining dynamics based on geometry" << std::endl;
    pcout << "Lattice " << lattice->getNx() << " by " << lattice->getNy() << " by " << lattice->getNz() << std::endl;
    pcout << "Parameters " << parameters->getNx() << " by " << parameters->getNy() << " by " << parameters->getNz() << std::endl;


    Box3D slice((plint) (parameters->getNx() / 2 + 1), (plint) (parameters->getNx() / 2 + 1),
                         0, parameters->getNy(),
                         0, parameters->getNz());

    // Loop over main time iteration.
    ValueTracer<T> converge(parameters->getLatticeU(), parameters->getResolution(), convergence_limit);

    pcout << "Timescale (iterations) = " << converge.getDeltaT() << std::endl;

    pcout << "Starting main loop" << std::endl;
    for (plint iT = 0; iT < parameters->nStep(maxT); ++iT)
    {
        converge.takeValue(getStoredAverageEnergy(*lattice), true);

        //Copied from ValueTracer.hasConverged
        //The original version checks for characteristic timescale,
        //which is too large (L / U / 2), about 3k iterations for a typical vessel
        //minIter is an arbitrary number.
        bool convergedFlag = 0;
        plint minIter = 20;
        if (iT > minIter && iT >= converge.getDeltaT()) {
            T average = converge.computeAverage();
            T stdDev = converge.computeStdDev(average);
            if (!std::isnan(stdDev/average))
            convergedFlag = (fabs(stdDev/average) < convergence_limit);
            else {
                logger() << "Simulation diverged." << std::endl;
                throw runtime_error("Simulation diverged.");
                break;
            }
        }

        if (iT % 500 == 0)
        {
            T maxspeed = computeMax(*computeVelocityComponent(*lattice, slice, 0), slice);
            pcout << "Max speed central slice: " << maxspeed << std::endl;
            computeStrainRateFromStress(*lattice, *wsstress, lattice->getBoundingBox());
            if(writeGifs)
                writeImage(*lattice, *parameters, *wsstress, iT, macroIter, imSize, slice);
        }

        if (convergedFlag)
        {
            pcout << "Flow solver has converged. This has lasted " << global::timer("main").stop() << " seconds." << endl;
            converge.resetValues();
            pcout << "Simulation converged for resolution = " << parameters->getResolution() << ", iterations = " << iT << std::endl;
            break;
        }

        global::timer("main").start();
        // Execute a time iteration.
        lattice->collideAndStream();
        global::timer("main").stop();
    }

    computeStrainRateFromStress(*lattice, *wsstress, lattice->getBoundingBox());
    if(writeGifs)
        writeImage(*lattice, *parameters, *wsstress, 999999, macroIter, imSize, slice);
}


void FlowSolver::execute(int argc, char* argv[])
{
    if (maskReceived)
    {
//                logger::info("geometry changed; recomputing flow");
        sourceToUnSerializer(new ReadFromSerialArray<int32_t> (geometry.data(), geometry.size()), charMask->getBlockUnSerializer(charMask->getBoundingBox(), IndexOrdering::forward));

        T uConst = 2.0; //default=2.0 for calculating max lattice speed. This value results in almost incompressible fluid
        generateParameters(uConst);
        printParameters(*lattice, *parameters);

        string path = instance->get_setting_as<std::string>("muscle_data_in_dir");
        string outPath = instance->get_setting_as<std::string>("muscle_data_out_dir");
        string flowFileName = instance->get_setting_as<std::string>("input_file");
        bool loadingSuccess = 0;
        if (macroIter == 0) {
//                    logger::info("setting up initial Poiseuille flow approximation...");
            initializeAtEquilibrium(*lattice, lattice->getBoundingBox(), CylindricalPoiseuilleDensityAndVelocity<T > (*parameters));

//                    logger::info("Attempting to load precalculated flow profile...");
            loadFlowFile(path, flowFileName, *lattice, loadingSuccess);

//                    logger::info("defining dynamics based on geometry");
            /// Boolmask defines isSolid. Fluid dynamics on fluid cells (isSolid = 0)

            plb::defineDynamics(*lattice, *charMask, new IncBGKdynamics<T, DESCRIPTOR>(parameters->getOmega()), 0);
            /// Boundary conditions are initialised here
            boundaryConditionSetup();
            /// No dynamics on all solid cells (isSolid = 1)
            plb::defineDynamics(*lattice, *charMask, new NoDynamics<T,DESCRIPTOR>, 1);
            //plb::defineDynamics(lattice, charMask, new BounceBack<T, DESCRIPTOR > (-3.), 1);
            /// Bounce back on all surface cells (isSolid = 2)
            plb::defineDynamics(*lattice, *charMask, new BounceBack<T, DESCRIPTOR > (-3.), 2);

//                    logger::info("Initializing LBM");
            lattice->initialize();

            if (!loadingSuccess) {
//                        logger::info("Using Poiseuille approximation for initial conditions.");
            }
        }
        else {
            /// No dynamics on all solid cells (isSolid = 1)
            plb::defineDynamics(*lattice, *charMask, new NoDynamics<T,DESCRIPTOR>, 1);
            /// Bounce back on all surface cells (isSolid = 2)
            plb::defineDynamics(*lattice, *charMask, new BounceBack<T, DESCRIPTOR > (-3.), 2);
            //don't override dynamics over BC cells
            plb::Box3D fluidDomain(1,parameters->getNx() - 2, 0, parameters->getNy() - 1, 0, parameters->getNz() - 1);
            plb::defineDynamics(*lattice, *charMask, fluidDomain, new IncBGKdynamics<T, DESCRIPTOR>(parameters->getOmega()), 0);
        }

        global::directories().setOutputDir(outPath);
        singleIteration();

        if (macroIter == 0 && !loadingSuccess) {
            global::directories().setOutputDir(path);
            saveBinaryBlock(*lattice, flowFileName);
        }

        bool writeVTKFile = true;
        if ((vtkIter == 0)||((macroIter % vtkIter) != 0))
            writeVTKFile = false;
        writeOutput(writeVTKFile);
    }
    else {
        // Send empty results back
        if (global::mpi().isMainProcessor()) {
            auto data = Data::nils(0);
            instance->send("shearStressOut", Message(timeInHours, data));
        }
        //logger::info("geometry not changed; sending same empty dataset");
    }
    macroIter++;
}

#endif
