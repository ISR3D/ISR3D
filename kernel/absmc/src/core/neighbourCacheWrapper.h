#ifndef ABSMC_NEIGHBOUR_CACHE_WRAPPER_H
#define ABSMC_NEIGHBOUR_CACHE_WRAPPER_H

#include <cstddef>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "core/neighbourCache.h"

namespace absmc {

/// Wrapper for NeighbourCache, offering convenience methods for neighbour searches.
template<class Agent>
class NeighbourCacheWrapper {
public:
    typedef Agent agent_t;

    /// For every agent contained in the agent vector, find all neighbours within range, and add pointers to them
    /// to the agent's neighbourhood record. If requested, the neighbourhood record is cleared before. Otherwise,
    /// duplicate entries are removed from the neighbourhood records afterwards.
    static void updateNeighbours(NeighbourCache<agent_t> const& nbCache, double range, bool doClear) {
        // implemented in terms of the method working on a single-agent
        std::vector<agent_t*> const& agents = nbCache.getAgentVector();
        const size_t nAgents = agents.size();
        #pragma omp parallel for schedule(dynamic, 100)
        for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
            updateNeighbours(iAgent, nbCache, range, doClear);
        }
    }

    /// Find all neighbours of the agent referenced by iAgent within range, and add pointers to them
    /// to the agent's neighbourhood record. If requested, the neighbourhood record is cleared before. Otherwise,
    /// duplicate entries are removed from the neighbourhood records afterwards.
    static void updateNeighbours(size_t iAgent, NeighbourCache<agent_t> const& nbCache, double range, bool doClear) {

        std::vector<agent_t*> const& agents = nbCache.getAgentVector();
        typename agent_t::nb_vec_t & nbPointers = agents[iAgent]->getNeighbours();

        // clear neighbourhood record if requested
        if (doClear) { nbPointers.clear(); }

        // find indices of all agents within a ball of radius range
        std::vector<size_t> nbIndices;
        nbCache.getNbIndices(iAgent, range, nbIndices);

        // add pointers to neighbourhood record
        for (size_t iIndex=0; iIndex<nbIndices.size(); iIndex++) {
            const size_t iAgent = nbIndices[iIndex];
            nbPointers.push_back(agents[iAgent]);
        }

        // remove possible duplicate entries from neighbourhood record
        // (neccessary only if record was not cleared before adding neigbours).
        // [could probably be implemented more efficient, but isn't expected to be executed frequently.]
        if (!doClear) {
            std::sort(nbPointers.begin(), nbPointers.end() );
            nbPointers.erase(std::unique(nbPointers.begin(), nbPointers.end()), nbPointers.end() );
        }
    }
};

} // end namespace absmc

#endif
