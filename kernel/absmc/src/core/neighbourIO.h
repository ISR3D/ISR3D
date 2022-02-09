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
    typedef AgentBase<nDim> agent_t;
    typedef typename AgentBase<nDim>::bond_vec_t bond_vec_t;

    static size_t writeFile(std::vector<index_vec_t> const& neighbourhoods, double length, std::string const& fileName);

    static size_t writeFile(AgentContainer<nDim> const& agentContainer, std::string const& fileName);

    // NBWriter does intentionally not support writing only a subset of agents contained in an AgentContainer.

private:
    static void findIndices();

    static void serializeNeighbours(std::vector<size_t> const& neighbours, std::vector<double> lengths, std::ostream & ostr);
};


template<size_t nDim>
class NBReader {
public:
    typedef std::vector<size_t> index_vec_t;
    typedef AgentBase<nDim> agent_t;
    typedef typename AgentBase<nDim>::bond_vec_t bond_vec_t;
    /// Read neighbourhood records from ascii text file.
    /// The user is responsible for matching the contents of the agent container with the correct file.
    static size_t readFile(std::string const& fileName, AgentContainer<nDim> & agentContainer, bool doClear=true);
};


} // end namespace absmc

#endif
