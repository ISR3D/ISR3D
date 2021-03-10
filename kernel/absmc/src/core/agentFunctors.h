#ifndef ABSMC_AGENT_FUNCTORS_H
#define ABSMC_AGENT_FUNCTORS_H

#include "core/geometry.h"
#include "core/agentRule.h"
#include "core/euclideanMetrics.h"

namespace absmc {

/// Functor setting the mobility vector of the second argument to the specified values.
template<class Agent>
class FSetMobility  : public std::binary_function<size_t, Agent, void> {
public:
    typedef Point<Agent::nDim, double> point_t;

    FSetMobility(point_t mobility_) : mobility (mobility_) { }
    void operator()(size_t id, Agent* agent) {
        agent->setMobility(mobility);
    }
private:
    point_t mobility;
};


/// Functor setting the initial position of the second argument to its current position.
template<class Agent>
class FResetPos0  : public std::binary_function<size_t, Agent, void> {
public:
    typedef Point<Agent::nDim, double> point_t;

    void operator()(size_t id, Agent* agent) {
        agent->setPos0(agent->getPos() );
    }
};


/// Functor setting the rule pointer of its second argument to the specified rule object,
/// assuming that the class passed as template parameter defines the setAgentRule(..) method.
/// Must be called with second argument of appropriate type only, e.g. via the typeId version of AgentContainer::forAll(..).
template<class Agent>
class FSetAgentRule : public std::binary_function<size_t, Agent, void> {
public:
    typedef Agent agent_t;
    typedef AgentBase<Agent::nDim> agent_base_t;
    typedef AgentRule<Agent> rule_t;

    FSetAgentRule(rule_t* rule_) : rule(rule_) { }

    void operator()(size_t id, agent_base_t* agentBase) {
        agent_t* agent = dynamic_cast<agent_t*>(agentBase);
        if (agent) { agent->setAgentRule(rule); }
    }
private:
    rule_t* rule;
};


/// Functor setting surface normals based on the number and location of neighbours close by.
/// Detection radius and max. number of neighbours to still be considered a surface cell are somewhat empiric.
template<class Agent>
class FDetectSurfaceFromNeighbours  : public std::binary_function<size_t, Agent, void> {
public:
    typedef Agent agent_t;
    typedef AgentBase<Agent::nDim> base_t;
    typedef Point<Agent::nDim, double> point_t;
    typedef EuclideanMetrics<agent_t::nDim, double> metrics_t;

    FDetectSurfaceFromNeighbours(double detectionRange_, size_t surfaceThreshold_,
                               NeighbourCache<base_t> ** nbCachePtr_)
    : detectionRange(detectionRange_), surfaceThreshold(surfaceThreshold_), nbCachePtr(nbCachePtr_) { }

    void operator()(size_t id, base_t * agent) {
        NeighbourCache<base_t>* nbCache = *nbCachePtr;
        assert(nbCache);

        std::vector<size_t> nbIndices;
        nbCache->getNbIndices(id, detectionRange, nbIndices);
        const size_t nNb = nbIndices.size();
        if (nNb > surfaceThreshold) {
            agent->unsetSurfaceNormal();
            return;
        }

        point_t accumulator;
        auto& agents = nbCache->getAgentVector();
        /// take average vector from nb to this point
        for (size_t iNb=0; iNb<nNb; iNb++) {
            accumulator -= agents[nbIndices[iNb] ]->getPos();
        }
        accumulator /= (double) nNb;
        accumulator += agent->getPos();

        /// normalize the vector
        accumulator /= norm<Agent::nDim>(accumulator);

        agent->setSurfaceNormal(accumulator);
    }

private:
    const double detectionRange;
    const size_t surfaceThreshold;
    NeighbourCache<base_t> ** nbCachePtr;
};

} // end namespace absmc

#endif
