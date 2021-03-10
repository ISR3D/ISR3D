#ifndef ABSMC_AGENT_FILE_READER_H
#define ABSMC_AGENT_FILE_READER_H

#include <cstddef>
#include <vector>
#include <util/agentTypeId.h>

#include "core/agentBase.h"
#include "core/agentFactory.h"
#include "core/agentContainer.h"

namespace absmc {

template<size_t nDim>
class AgentFileReader {
public:
    typedef AgentBase<nDim> agent_t;
    typedef Point<nDim, double> point_t;
    typedef std::vector<agent_t*> nb_vec_t;

    /// Read agent data from ascii file, create agents using the given factory, and add them to the given vector.
    /// If typeId is specified, agents are created only for the matching records.
    /// File format: x,y,x0,y0,r separated by whitespace.
    static size_t readFile(std::string const& fileName,
                           AgentFactory<nDim> const& factory,
                           std::vector<agent_t*> & agents,
                           AgentTypeId typeId=tAny);

    /// Read agent data from ascii file, create agents using the given factory, and add them to the given agent container.
    /// If typeId is specified, agents are created only for the matching records.
    /// File format: x,y,x0,y0,r separated by whitespace.
    static size_t readFile(std::string const& fileName,
                           AgentFactory<nDim> const& factory,
                           AgentContainer<nDim> & agentContainer,
                           AgentTypeId typeId=tAny);

    /// Read initial agent data from ascii csv file, create agents using the given factory, and add them to the given vector.
    /// File format: x,y,z,r separated by whitespace.
    static size_t readCsvFile(std::string const& fileName,
                          AgentFactory<nDim> const& factory,
                          std::vector<agent_t*> & agents);


    /// Read initial agent data from ascii csv file, create agents using the given factory, and add them to the given agent container.
    /// File format: x,y,z,r separated by whitespace.
    static size_t readCsvFile(std::string const& fileName,
                          AgentFactory<nDim> const& factory,
                          AgentContainer<nDim> & agentContainer);


};

} // end namespace absmc

#endif
