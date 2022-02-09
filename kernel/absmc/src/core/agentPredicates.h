#ifndef ABSMC_AGENT_PREDICATES_H
#define ABSMC_AGENT_PREDICATES_H

#include <util/agentTypeId.h>
#include "core/geometry.h"
#include <random>

namespace absmc {

/// Predicate returning true if agent coordinate iDim is closer than range to the [xLo,xHi] interval boundary or outside the interval.
template<class Agent>
class PSelectCloseToBoundary : public std::binary_function<size_t, Agent, bool> {
public:
    PSelectCloseToBoundary(size_t iDim_, double xLo_, double xHi_, double range_)
        : iDim(iDim_), xLo(xLo_), xHi(xHi_), range(range_) { }
    bool operator()(size_t id, Agent* agent) {
        const double x = agent->getPos()[iDim];
        return ((x-xLo < range) || (xHi-x < range));
    }
private:
    size_t iDim;
    double xLo, xHi, range;
};

/// Predicate returning true if agent coordinate iDim is further than range to the [xLo,xHi] interval boundary.
template<class Agent>
class PSelectAwayFromBoundary : public std::binary_function<size_t, Agent, bool> {
public:
    PSelectAwayFromBoundary(size_t iDim_, double xLo_, double xHi_, double range_)
        : iDim(iDim_), xLo(xLo_), xHi(xHi_), range(range_) { }
    bool operator()(size_t id, Agent* agent) {
        const double x = agent->getPos()[iDim];
        return ((x-xLo > range) && (xHi-x > range));
    }
private:
    size_t iDim;
    double xLo, xHi, range;
};

/// Predicate returning true if agent is inside [xLo,xHi] interval.
template<class Agent>
class PSelectInterval : public std::binary_function<size_t, Agent, bool> {
public:
    PSelectInterval(size_t iDim_, double xLo_, double xHi_)
        : iDim(iDim_), xLo(xLo_), xHi(xHi_) { }
    bool operator()(size_t id, Agent* agent) {
        const double x = agent->getPos()[iDim];
        return ((x >= xLo) && (xHi >= x));
    }
private:
    size_t iDim;
    double xLo, xHi;
};

/// Predicate returning true if agent is inside [xLo,xHi] interval.
template<class Agent>
class PSelectOutsideInterval : public std::binary_function<size_t, Agent, bool> {
public:
    PSelectOutsideInterval(size_t iDim_, double xLo_, double xHi_)
        : iDim(iDim_), xLo(xLo_), xHi(xHi_) { }
    bool operator()(size_t id, Agent* agent) {
        const double x = agent->getPos()[iDim];
        return ((x <= xLo) || (xHi <= x));
    }
private:
    size_t iDim;
    double xLo, xHi;
};

/// Predicate returning true if agent is more than L away from pos0.
template<class Agent>
class PSelectByDisplacement : public std::binary_function<size_t, Agent, bool> {
public:
    PSelectByDisplacement(double L_)
        : L(L_) { }
    bool operator()(size_t id, Agent* agent) {
        const Point<3, double> pos = agent->getPos();
        const Point<3, double> pos0 = agent->getPos0();
        return (norm<3>(pos0 - pos) > L);
    }
private:
    double L;
};

/// Predicate returning true for the specified fraction of randomly chosen agents.
template<class Agent>
class PSelectRandomly : public std::binary_function<size_t, Agent, bool> {
public:
    PSelectRandomly(double fractionToSelect_) : fractionToSelect(fractionToSelect_), generator(), distribution(0., 1.){ }

    bool operator()(size_t id, Agent* agent) {
        return (distribution(generator) < fractionToSelect);
    }
private:
    double fractionToSelect;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;
};


/// Predicate returning true for agents of specified type
template<class Agent>
class PSelectAgentsByType : public std::binary_function<size_t, Agent, bool> {
public:
    typedef Point<Agent::nDim, double> point_t;

    PSelectAgentsByType( AgentTypeId typeId_)
        :  typeId(typeId_) { }

    bool operator()(size_t id, Agent* agent) {
        // filter out agents not matching the specified typeId
        if (typeId != tAny && agent->getTypeId() != typeId) return false;
        return true;
    }
private:
    const AgentTypeId typeId;
};

} // end namespace absmc


#endif
