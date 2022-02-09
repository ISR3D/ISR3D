/// This program loads an artery from a .csv file,
/// sets up types and bonds and
/// outputs a .dat file suitable for equilibration (stage1)

#include <iostream>

#include "absmc3D.h"

using namespace absmc;

const std::string solverName = "csvArteryLoader3D";

/// main
int main(int argc, char *argv[])
{
    if (argc!=2) {
        std::cout << "Usage: " << solverName << " configDir/configName" << std::endl;
        return EXIT_FAILURE;
    }

    const std::string configFileName(argv[1]);

    AgentSolver<3> slv(solverName);
    slv.loadConfig(configFileName);
    slv.readIOFileNames("csv_input_file", "csv_output_file");
    slv.readCharLengthScales();

    // load agents from file and put them in AgentContainer
    slv.readCsv3DFile(slv.getInFileName());

    // create neighbour detector and initialise neighbourhoods
    slv.initStaticNeighbourDetector();

    // create bonds from the detected neighbourhoods
    slv.makeBondsFromNeighbours();

    // write the final configuration to file
    slv.readVtpParameters("reader_vtp_scalars", "reader_vtp_vectors");
    slv.writeAllOutputs("stage1.");

    return EXIT_SUCCESS;
}
