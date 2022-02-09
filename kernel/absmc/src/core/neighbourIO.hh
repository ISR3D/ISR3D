#ifndef ABSMC_NEIGHBOUR_IO_HH
#define ABSMC_NEIGHBOUR_IO_HH

#include <algorithm>
#include <cassert>

#include "core/neighbourIO.h"
#include "core/agentBase.h"
#include <util/miscUtil.h>

namespace absmc {

template<size_t nDim>
size_t NBWriter<nDim>::writeFile(std::vector<index_vec_t> const& neighbourhoods, double length, std::string const& fileName)
{
    std::ofstream ofstr (fileName.c_str());
    if (!ofstr) {
        std::cout << "Cannot open file " << fileName << " for writing." << std::endl;
        return 0;
    }

    const size_t numNbh = neighbourhoods.size();
    for (size_t iNbh=0; iNbh<numNbh; iNbh++) {
        index_vec_t const& curNeighbourhood = neighbourhoods[iNbh];

        const size_t numNb= curNeighbourhood.size();
        std::vector<double> lengths;
        lengths.assign (numNb,length);

        serializeNeighbours(curNeighbourhood, lengths, ofstr);
    }
    return numNbh;
}

template<size_t nDim>
size_t NBWriter<nDim>::writeFile(AgentContainer<nDim> const& agentContainer, std::string const& fileName)
{
    std::ofstream ofstr (fileName.c_str());
    if (!ofstr) {
        std::cout << "Cannot open file " << fileName << " for writing." << std::endl;
        return 0;
    }

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

        bond_vec_t const& nbPointers = agents[iAgent]->getBonds();
        index_vec_t nbIndices;
        std::vector<double> lengths;

        const size_t numNb= nbPointers.size();
        for (size_t iNb=0; iNb<numNb; iNb++) {
            agent_t* curP = nbPointers[iNb].other;
            const size_t curId = curP->getId();
            nbIndices.push_back(agentPos[curId]);
            lengths.push_back(nbPointers[iNb].length0);
        }

        serializeNeighbours(nbIndices, lengths, ofstr);
    }
    return nAgents;
}

template<size_t nDim>
void NBWriter<nDim>::serializeNeighbours(index_vec_t const& neighbours, std::vector<double> lengths, std::ostream & ostr)
{
    const size_t numNb= neighbours.size();
    /// write indices, separated by colon
    for (size_t iNb=0; iNb<numNb; iNb++) {
        ostr  << neighbours[iNb] << ":";
    }
    /// always write a least one colon, so that there are no empty lines
    if (numNb==0) { ostr << ":"; }
    ostr << std::endl;

    /// same for bond lengths
    for (size_t iNb=0; iNb<numNb; iNb++) {
        ostr  << std::setprecision (5) << lengths[iNb] << ":";
    }
    /// always write a least one colon, so that there are no empty lines
    if (numNb==0) { ostr << ":"; }
    ostr << std::endl;
}

template<size_t nDim>
size_t NBReader<nDim>::readFile(std::string const& fileName, AgentContainer<nDim> & agentContainer, bool doClear)
{
    std::ifstream ifstr(fileName.c_str());
    if (!ifstr) {
        std::cout << "Cannot open file " << fileName << " for reading."  << std::endl;
        return 0;
    }

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

        std::getline(ifstr, line);
        std::vector<double> lengths;
        util::tokenizeAndConvert<double>(line, lengths, ':');

        /// get reference to current agent's neighbourhood
        bond_vec_t & nbPointers = agents[numRead]->getBonds();
        if (doClear) { nbPointers.clear(); }

        /// add pointers to neighbourhood of current agent
        for (size_t iNb=0; iNb<numNb; iNb++) {
            agents[numRead]->addBond(agents[nbIndices[iNb]], lengths[iNb], false);
        }
        numRead++;
    }

    return numRead;
}


} // end namespace absmc

#endif
