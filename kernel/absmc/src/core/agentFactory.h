#ifndef ABSMC_AGENT_FACTORY_H
#define ABSMC_AGENT_FACTORY_H

#include <cstddef>
#include <vector>

#include <util/agentTypeId.h>

#include "core/agentBase.h"
#include "core/geometry.h"
#include "core/agentContainer.h"

namespace absmc {

/// Base class for agent factories.
template<size_t nDim>
class AgentFactory {
public:
    typedef AgentBase<nDim> agent_t;
    typedef Point<nDim, double> point_t;

    virtual ~AgentFactory() { }

    virtual agent_t* create(AgentTypeId typeId, point_t pos, point_t pos0, double r) const = 0;
    virtual agent_t* create(AgentTypeId typeId, point_t pos, point_t pos0, double r,
                            double drugConc,  double wssOsi, double wssMax, double neoDist, double nitricConc,
                            int age, int clock, int lengthG1,
                            bool baFlag, bool ciFlag, bool selectedForNOFlag) const = 0;
};


/// Create agents from the given vector of points using the given factory, assuming constant type and radius.
template<int nDim>
void pointsToCells(std::vector<Point<nDim, double> > const& points,
                   AgentContainer<nDim> & agentContainer,
                   AgentFactory<nDim> const& agentFactory, AgentTypeId typeId, double r)
{
    const size_t nPoints = points.size();
    for (size_t iPoint=0; iPoint<nPoints; iPoint++) {
        AgentBase<nDim>* newAgent = agentFactory.create(typeId, points[iPoint], points[iPoint], r);
        if (newAgent) { agentContainer.add(newAgent); }
    }
}

} // end namespace absmc

#endif
