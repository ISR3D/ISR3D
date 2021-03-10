#ifndef ABSMC_CURRENT_NEIGHBOUR_DETECTOR_H
#define ABSMC_CURRENT_NEIGHBOUR_DETECTOR_H

#include <cstddef>
#include <cassert>
#include <vector>

#include "core/neighbourDetector.h"
#include "core/neighbourCache.h"
#include "core/neighbourCacheWrapper.h"

namespace absmc {

/// "Current" neighbourhood detection policy: compute neighbourhood records from scratch at every call to update().
/// All agents within a sphere of radius range at the time of updating are considered neighbours.
/// Attention, this is a rather expensive policy! You will most probably want to use the VerletNeighbourDetector
/// in conjunction with a force using a cut-off rang instead.
template<class Agent>
class CurrentNeighbourDetector : public NeighbourDetector<Agent> {
public:
    typedef Agent agent_t;

    CurrentNeighbourDetector(double range_) : range(range_) { }

    /// Initialise neighbours for all agents contained in the agent vector.
    virtual void init(std::vector<agent_t*> agents) {
        NeighbourCache<agent_t> nbCache(agents);
        NeighbourCacheWrapper<agent_t>::updateNeighbours(nbCache, range, true);
    }

    /// Update neighbours for all agent_t contained in the agent vector,
    /// on a slow time scale. Called after deployment step and after application of biological ruleset.
    virtual void updateSlowTimeScale(std::vector<agent_t*> agents) {
        NeighbourCache<agent_t> nbCache(agents);
        NeighbourCacheWrapper<agent_t>::updateNeighbours(nbCache, range, true);
    }

    /// Update neighbours for all agent_t contained in the agent vector,
    /// on a fast time scale. Called at every iteration of the integrator.
    virtual void updateFastTimeScale(std::vector<agent_t*> agents, double displMax) {
        NeighbourCache<agent_t> nbCache(agents);
        NeighbourCacheWrapper<agent_t>::updateNeighbours(nbCache, range, true);
    }

private:
    const double range;   //  interaction range
};

} // end namespace absmc

#endif
