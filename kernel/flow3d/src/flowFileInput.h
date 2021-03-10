#ifndef FLOWFILEINPUT_H
#define FLOWFILEINPUT_H

#include <iostream>
#include <fstream>
#include <string>

void loadFlowFile(std::string path, std::string flowFileName,
        MultiBlockLattice3D<T, DESCRIPTOR>& lattice, bool& loadingSuccess)
{

    global::directories().setOutputDir(path);
    try
    {
        loadBinaryBlock(lattice, flowFileName);
        loadingSuccess = 1;
        pcout << "Precalculated flow profile loaded." << std::endl;
    } catch (...)
    {
        loadingSuccess = 0;
        pcout << "Failed to load precalculated flow." << std::endl;
    }

    global::mpi().bCast(&loadingSuccess, 1);
}

#endif // FLOWFILEINPUT_H
