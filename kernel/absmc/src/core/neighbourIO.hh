#ifndef ABSMC_NEIGHBOUR_IO_HH
#define ABSMC_NEIGHBOUR_IO_HH

#include <algorithm>
#include <cassert>

#include "core/neighbourIO.h"
#include "core/agentBase.h"
#include <util/miscUtil.h>

namespace absmc {

template<size_t nDim>
size_t NBWriter<nDim>::writeFile(std::vector<index_vec_t> const& neighbourhoods, std::string const& fileName)
{
    std::ofstream ofstr (fileName.c_str());
    if (!ofstr) {
        std::cout << "Cannot open file " << fileName << " for writing." << std::endl;
        return 0;
    }

    const size_t numNbh = neighbourhoods.size();
    for (size_t iNbh=0; iNbh<numNbh; iNbh++) {
        index_vec_t const& curNeighbourhood = neighbourhoods[iNbh];
        serializeNeighbours(curNeighbourhood, ofstr);
    }
    return numNbh;
}

template<size_t nDim>
size_t NBWriter<nDim>::writeFile(AgentContainer<nDim> const& agentContainer, std::string const& fileName, bool writeBonds)
{
    std::ofstream ofstr (fileName.c_str());
    if (!ofstr) {
        std::cout << "Cannot open file " << fileName << " for writing." << std::endl;
        return 0;
    }

    typedef AgentBase<nDim> agent_t;
    std::vector<agent_t*> const& agents = agentContainer.getAgentVector();
    const size_t nAgents = agents.size();

    /// this vector stores agent position in the container for each agentId.
    /// procedure becomes ineffective if there are lots of deletions in the simulation
    std::vector<size_t> agentPos;
    agentPos.resize(2 * nAgents);

    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        const size_t id = agents[iAgent]->getId();
        if (id > agentPos.size())
            agentPos.resize(2 * id);
        agentPos[id] = iAgent;
    }

    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {

        std::vector<agent_t*> const& nbPointers = writeBonds ? agents[iAgent]->getBonds() : agents[iAgent]->getNeighbours();
        index_vec_t nbIndices;

        const size_t numNb= nbPointers.size();
        for (size_t iNb=0; iNb<numNb; iNb++) {
            agent_t* curP = nbPointers[iNb];
            const size_t curId = curP->getId();
            nbIndices.push_back(agentPos[curId]);
        }

        serializeNeighbours(nbIndices, ofstr);
    }
    return nAgents;
}

template<size_t nDim>
void NBWriter<nDim>::serializeNeighbours(index_vec_t const& neighbours, std::ostream & ostr)
{
    const size_t numNb= neighbours.size();
    /// write indices, separated by colon
    for (size_t iNb=0; iNb<numNb; iNb++) {
        ostr  << neighbours[iNb] << ":";
    }
    /// always write a least one colon, so that there are no empty lines
    if (numNb==0) { ostr << ":"; }
    ostr << std::endl;
}

template<size_t nDim>
size_t NBReader<nDim>::readFile(std::string const& fileName, AgentContainer<nDim> & agentContainer, bool doClear, bool readBonds)
{
    std::ifstream ifstr(fileName.c_str());
    if (!ifstr) {
        std::cout << "Cannot open file " << fileName << " for reading."  << std::endl;
        return 0;
    }

    typedef AgentBase<nDim> agent_t;
    typedef typename AgentBase<nDim>::nb_vec_t nb_vec_t;
    typedef std::vector<size_t> index_vec_t;

    std::vector<agent_t*> const& agents = agentContainer.getAgentVector();

    // TODO: implement simple checks

    size_t numRead = 0;
    while (ifstr) {
        std::string line;
        std::getline(ifstr, line);

        /// an empty line indicates end of file,
        /// since in every line corresponding to an agent there must be at least one colon.
        if (line.empty() ) break;

        /// unserialize one neighbourhood record
        index_vec_t nbIndices;
        util::tokenizeAndConvert<size_t>(line, nbIndices, ':');
        const size_t numNb = nbIndices.size();

        /// get reference to current agent's neighbourhood
        assert(numRead < agents.size());
        nb_vec_t & nbPointers = readBonds ? agents[numRead]->getBonds() : agents[numRead]->getNeighbours();
        if (doClear) { nbPointers.clear(); }

        /// add pointers to neighbourhood of current agent
        for (size_t iNb=0; iNb<numNb; iNb++) {
            assert(nbIndices[iNb] < agents.size());
            nbPointers.push_back(agents[nbIndices[iNb]] );
        }
        numRead++;
    }

    return numRead;
}


} // end namespace absmc

#endif
