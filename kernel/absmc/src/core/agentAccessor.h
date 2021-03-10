#ifndef ABSMC_AGENT_ACCESSOR_H
#define ABSMC_AGENT_ACCESSOR_H

#include <vector>

#include "core/geometry.h"

namespace absmc {

/// Accessor class used by kd-tree to access spatial coordinates of agents.
template<class Agent>
class AgentAccessor {
public:
    typedef Agent agent_t;
    typedef Point<Agent::nDim, double> point_t;

    /// constructor
    AgentAccessor(std::vector<agent_t*> const& data_) : data(data_) { }
    // compiler-generated copy constructor and assignment operator are fine.

    /// return number of objects stored in the underlying container
    size_t size() const {
        return data.size();
    }

    /// get coordinate iDim of element iData.
    double getCoordinate(size_t iData, size_t iDim) const {
        assert(iData < size());
        assert(iDim < Agent::nDim);
        return (data[iData]->getPos())[iDim];
    }

    /// set coordinate iDim of element iData.
    void setCoordinate(size_t iData, size_t iDim, double value) {
        assert(iData < size());
        assert(iDim < Agent::nDim);
        (data[iData]->getPos())[iDim] = value;
    }

    /// construct a Point data structure from spatial data accessible by acc.
    point_t makePoint(size_t iData) const {
        assert(iData < size());
        return data[iData]->getPos();
    }

private:
    std::vector<agent_t*> const& data;
};

} // end namespace absmc

#endif
