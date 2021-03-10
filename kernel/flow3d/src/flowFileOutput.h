#ifndef FLOWFILEOUTPUT_H
#define FLOWFILEOUTPUT_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

template<typename T>
std::string toString(const T& value)
{
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

void printParameters(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, IncomprFlowParam<T> const& parameters)
{
    pcout << getMultiBlockInfo(lattice);
    pcout << "Delta T=" << parameters.getDeltaT() << endl;
    pcout << "Delta X=" << parameters.getDeltaX() << endl;
    pcout << "Viscosity in lattice units=" << parameters.getLatticeNu() << endl;
    pcout << "Velocity in lattice units (proportional to Mach number) =" << parameters.getLatticeU() << endl;
    pcout << "x Length =" << parameters.getLx() << endl;
    pcout << "y Length =" << parameters.getLy() << endl;
    pcout << "z Length =" << parameters.getLz() << endl;
    pcout << "Omega =" << parameters.getOmega() << endl;
    pcout << "Reynolds =" << parameters.getRe() << endl;
    pcout << "Resolution =" << parameters.getResolution() << endl;
    pcout << "Relaxation Time =" << parameters.getTau() << endl;
    pcout << "conversion from dimensionless to lattice units for space coordinate  =" << parameters.nCell((T) 1) << endl;
    pcout << "conversion from dimensionless to lattice units for time coordinate   =" << parameters.nStep((T) 1) << endl;
    pcout << "Parameters Nx =" << parameters.getNx() << endl;
    pcout << "Parameters Ny =" << parameters.getNy() << endl;
    pcout << "Parameters Nz =" << parameters.getNz() << endl;
    pcout << "Lattice Nx =" << lattice.getNx() << endl;
    pcout << "Lattice Ny =" << lattice.getNy() << endl;
    pcout << "Lattice Nz =" << lattice.getNz() << endl;
}

void writeVTK(MultiBlockLattice3D<T, DESCRIPTOR>& lattice,
              IncomprFlowParam<double> const& parameters, plint iter)
{
    double dx = parameters.getDeltaX();
    double dt = parameters.getDeltaT();
    VtkImageOutput3D<double> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx / dt);
    vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.);
    //vtkOut.writeData<float>(*subtract(*computeDensity(lattice),(T)1.), "density", 1.);
    vtkOut.writeData < 3, float>(*computeVelocity(lattice), "velocity", dx / dt);
}

void writeImage(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, IncomprFlowParam<double> const& parameters,
                MultiTensorField3D<T, 6 > & wsstress,
                //MultiScalarField3D<T > & wss,
                plint iT, plint bigIter, plint imSize, Box3D & slice)
{
    T dt = parameters.getDeltaT();
    //write .gifs (YZ plane)
    pcout << "Writing image at dimensionless time " << iT * dt << std::endl;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(
                               createFileName("velocity"+toString(bigIter)+"_", iT, 6),
                               *computeVelocityNorm(lattice, slice), imSize, imSize);
    imageWriter.writeScaledGif(
                               createFileName("density"+toString(bigIter)+"_", iT, 6),
                               *computeDensity(lattice, slice), imSize, imSize);
    imageWriter.writeScaledGif(
                               createFileName("wsstress"+toString(bigIter)+"_", iT, 6),
                               *multiply(10. /* * shearStressConversionFactor(parameters) */, *computeSymmetricTensorNorm(wsstress, slice), lattice.getBoundingBox())
                               , imSize, imSize);
//    imageWriter.writeScaledGif(
//            createFileName("wss", iT + bigIT, 6),wss, imSize, imSize);

    const plint Nx = parameters.getNx();
    const plint Ny = parameters.getNy();
    const plint Nz = parameters.getNz();
    //logger::info("Nx %d, Ny %d, Nz %d", Nx, Ny, Nz);

    Box3D slicey(0, Nx, (plint) (Ny / 2 + 1), (plint) (Ny / 2 + 1), 0, Nz);
    Box3D slicez(0, Nx, 0, Ny, (plint) (Nz / 2 + 1), (plint) (Nz / 2 + 1));

    Box3D zaxis((plint) (Nx / 2 + 1), (plint) (Nx / 2 + 1), (plint) (Ny / 2 + 1), (plint) (Ny / 2 + 1), 0, Nz);
    Box3D yaxis((plint) (Nx / 2 + 1), (plint) (Nx / 2 + 1), 0, Ny, (plint) (Nz / 2 + 1), (plint) (Nz / 2 + 1));


    //write .gifs (XZ and XY planes)
    imageWriter.writeScaledGif(
                               createFileName("velocity_y"+toString(bigIter)+"_", iT, 6),
                               *computeVelocityNorm(lattice, slicey), Nx*2, Nz*2);
    imageWriter.writeScaledGif(
                               createFileName("velocity_z"+toString(bigIter)+"_", iT, 6),
                               *computeVelocityNorm(lattice, slicez), Nx*2, Ny*2);

    imageWriter.writeScaledGif(
                               createFileName("density_y"+toString(bigIter)+"_", iT, 6),
                               *computeDensity(lattice, slicey), Nx*2, Ny*2);


    //write velocity values along Y and Z axes on the middle slice
    MultiScalarField3D<T> velocity_yaxis = *computeVelocityNorm(lattice, yaxis);
    plb_ofstream ofilevel_y("vel_yaxis_mid.dat");
    ofilevel_y << velocity_yaxis << endl;
    MultiScalarField3D<T> velocity_zaxis = *computeVelocityNorm(lattice, zaxis);
    plb_ofstream ofilevel_z("vel_zaxis_mid.dat");
    ofilevel_z << velocity_zaxis << endl;
}

#endif // FLOWFILEOUTPUT_H
