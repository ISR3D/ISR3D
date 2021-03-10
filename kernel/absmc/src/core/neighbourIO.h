#ifndef ABSMC_NEIGHBOUR_IO_H
#define ABSMC_NEIGHBOUR_IO_H

#include <cstddef>
#include <vector>
#include <fstream>

#include <util/agentTypeId.h>
#include "core/agentContainer.h"

namespace absmc {

template<size_t nDim>
class NBWriter {
public:
    typedef std::vector<size_t> index_vec_t;

    static size_t writeFile(std::vector<index_vec_t> const& neighbourhoods, std::string const& fileName);

    static size_t writeFile(AgentContainer<nDim> const& agentContainer, std::string const& fileName, bool writeBonds = false);

    // NBWriter does intentionally not support writing only a subset of agents contained in an AgentContainer.

private:
    static void findIndices();

    static void serializeNeighbours(std::vector<size_t> const& neighbours, std::ostream & ostr);
};


template<size_t nDim>
class NBReader {
public:
    /// Read neighbourhood records from ascii text file.
    /// The user is responsible for matching the contents of the agent container with the correct file.
    static size_t readFile(std::string const& fileName, AgentContainer<nDim> & agentContainer, bool doClear=true, bool readBonds = false);
};


} // end namespace absmc

#endif
