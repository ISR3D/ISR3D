#ifndef ABSMC_NEIGHBOUR_CACHE_H
#define ABSMC_NEIGHBOUR_CACHE_H

#include <cstddef>
#include <vector>
#include <algorithm>

#include "kdtree/kdtree.h"
#include "core/agentAccessor.h"

namespace absmc {

/// The NeighbourCache class offers a low-level method for finding neighbours of an agent within the given agent vector.
/// To make this method fast, it builds a structure for fast spatial searches at construction; the actual implementation
/// encapsulates a kd-tree.
template<class Agent>
class NeighbourCache {
public:
    typedef Agent agent_t;

    /// Build cache for fast neighbour searches; once the cache has been built, any change of position of the agents
    /// contained in the underlying agent vector invalidate the cache. The user is responsible to take care of this.
    NeighbourCache(std::vector<agent_t*> const& agents_) : agents(agents_), acc(agents), kdtree(acc) { }

    /// Get a reference to the vector of agents that the cache is constructed on top of.
    std::vector<agent_t*> const& getAgentVector() const { return agents; }

    /// Low-level method to retrieve indices (into the agent vector) of all neighbours of the agent
    /// referenced by iAgent. A neighbour is defined as another agent within a ball of radius range.
    void getNbIndices(size_t iAgent, double range, std::vector<size_t> & nbIndices) const {
        nbIndices.clear();
        kdtree.findWithinBall(iAgent, range, nbIndices);
    }

private:
    std::vector<agent_t*> const& agents;
    const AgentAccessor<agent_t> acc;
    const kdtree::KDTree<Agent::nDim, double, AgentAccessor<agent_t> > kdtree;
};

} // end namespace absmc

#endif
