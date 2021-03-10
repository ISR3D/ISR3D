#ifndef ABSMC_STATIC_NEIGHBOUR_DETECTOR_H
#define ABSMC_STATIC_NEIGHBOUR_DETECTOR_H

#include <cstddef>
#include <vector>

#include "core/neighbourDetector.h"
#include "core/neighbourCache.h"
#include "core/neighbourCacheWrapper.h"

namespace absmc {

/// "Static" neighbourhood detection policy: neighbourhoods are determined at initialisation.
/// All agents within a sphere of radius range at the time of initialisation are considered neighbours.
template<class Agent>
class StaticNeighbourDetector : public NeighbourDetector<Agent> {
public:
    typedef Agent agent_t;

    StaticNeighbourDetector(double range_) : range(range_) { }

    /// Initialise neighbours for all agents contained in the agent vector;
    virtual void init(std::vector<agent_t*> agents) {
        NeighbourCache<agent_t> nbCache(agents);
        NeighbourCacheWrapper<agent_t>::updateNeighbours(nbCache, range, true);
    }

    /// Update neighbours for all agent_t contained in the agent vector,
    /// on a slow time scale. Called after deployment step and after application of biological ruleset.
    virtual void updateSlowTimeScale(std::vector<agent_t*> agents) {
        // nothing to do, since neighbourhood records are static.
    }

    /// Update neighbours for all agent_t contained in the agent vector,
    /// on a fast time scale. Called at every iteration of the integrator.
    virtual void updateFastTimeScale(std::vector<agent_t*> agents, double displMax) {
        // nothing to do, since neighbourhood records are static.
    }

private:
    double range;

};


} // end namespace absmc

#endif
