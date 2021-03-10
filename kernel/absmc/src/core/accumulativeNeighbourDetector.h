#ifndef ABSMC_ACCUMULATIVE_NEIGHBOUR_DETECTOR_H
#define ABSMC_ACCUMULATIVE_NEIGHBOUR_DETECTOR_H

#include <cstddef>
#include <vector>

#include "core/neighbourDetector.h"
#include "core/neighbourCache.h"
#include "core/neighbourCacheWrapper.h"

namespace absmc {

/// "Accumulative" neighbourhood detection policy: at initialisation, all agents within a sphere of radius range
/// are considered neighbours. At slow time scale updates, existing neighbours are kept, and agents newly detected
/// within range are added. Neighbourhood records are not changed at fast time scale updates.
template<class Agent>
class AccumulativeNeighbourDetector : public NeighbourDetector<Agent> {
public:
    typedef Agent agent_t;

    AccumulativeNeighbourDetector(double range_) : range(range_) { }

    /// Initialise neighbours for all agents contained in the agent vector.
    virtual void init(std::vector<agent_t*> agents) {
        NeighbourCache<agent_t> nbCache(agents);
        NeighbourCacheWrapper<agent_t>::updateNeighbours(nbCache, range, true);
    }

     /// Update neighbours for all agent_t contained in the agent vector,
    /// on a slow time scale. Called after deployment step and after application of biological ruleset.
    virtual void updateSlowTimeScale(std::vector<agent_t*> agents) {
        // add new neighbours within range, without dropping those who might have moved out of range.
        NeighbourCache<agent_t> nbCache(agents);
        NeighbourCacheWrapper<agent_t>::updateNeighbours(nbCache, range, false);
    }

    /// Update neighbours for all agent_t contained in the agent vector,
    /// on a fast time scale. Called at every iteration of the integrator.
    virtual void updateFastTimeScale(std::vector<agent_t*> agents, double displMax) {
        // no neighbour accumulation on fast time scale: nothing to do.
    }

private:
    const double range;
};


} // end namespace absmc

#endif
