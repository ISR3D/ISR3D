/// This program loads an artery from a .csv file,
/// sets up types and bonds and
/// outputs a .dat file suitable for equilibration (stage1)

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

using std::cout;
using std::cerr;
using std::endl;

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
    const std::string outputDir = createFileName(configDir, "data.csvreader", "");

    config().readFile(configFileName);

    // agent input file
    const std::string inFileName = config().getValue<std::string>("csv_input_file");
    const std::string outFileName = config().getValue<std::string>("csv_output_file");

    const auto fileNameParts = parseFileName(outFileName);
    const std::string baseName = fileNameParts.baseName;

    // characteristic length scale
    const double L = config().getValue<double>("smc_mean_rad");
    const double sigma = config().getValue<double>("smc_rad_sigma");

    // load agents from file includes smc+iel+obstacle agents
    AgentContainer<3> agents;
    AgentFileReader<3>::readCsvFile(createFileName(configDir, inFileName), AgentFactory3D(), agents);
    cout << agents.count() << " agents loaded from csv." << endl;

    // create neighbour detector and initialise neighbourhoods
    StaticNeighbourDetector<agent_t> nbDetector(2.2 * L + 2 * sigma);
    nbDetector.init(agents.getAgentVector());

    const size_t nAgents = agents.count();

    // create bonds from the detected neighbourhoods
    /// a bit dirty, can be a container-wide thing or pre-set in csv
    for (size_t iAgent = 0; iAgent < nAgents; iAgent++) {
        agents.getAgentVector()[iAgent]->addNeighboursToBonds();
    }
    cout << "Bonds initialized from neighbours" << endl;

    // NeighbourCache<agent_t>* nbCacheForCiRule = new NeighbourCache<agent_t > (agents.getAgentVector());
    // DetectSurfaceFromNeighboursRule3D surfaceDetectorRule(L * 11, 100, &nbCacheForCiRule);
    // agents.forAll(FSetAgentRule<CellBase3D>(&surfaceDetectorRule) );
    // agents.execAgentRules();

    // vtp output parameters
    const std::string vtpScalars = config().getValue<std::string>("reader_vtp_scalars");
    const std::string vtpVectors = config().getValue<std::string>("reader_vtp_vectors");

    // write the final configuration to file
    agents.writeFile(createFileName(configDir, "stage1."+baseName, "dat") );
    NBWriter<3>::writeFile(agents, createFileName(configDir, "stage1."+baseName+"_nb", "dat") );
    VtkWriter<agent_t>::writeVtkPolyData(agents.getAgentVector(), createFileName(configDir, "stage1."+baseName, "vtp"), vtpScalars, vtpVectors);

    return EXIT_SUCCESS;
}
