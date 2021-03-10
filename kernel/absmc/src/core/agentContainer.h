#ifndef ABSMC_AGENT_CONTAINER_H
#define ABSMC_AGENT_CONTAINER_H

#include <cstddef>
#include <vector>
#include <set>

#include <util/agentTypeId.h>
#include "core/agentBase.h"

namespace absmc {

template<size_t nDim>
class AgentContainer {
public:

    typedef AgentBase<nDim> agent_t;
    typedef Point<3, double> point_t;

    /// Constructor: reserve memory for pointers.
    AgentContainer(size_t sizeToReserve=4096) : agents() { agents.reserve(sizeToReserve); }

    /// Destructor: delete all agents contained.
    ~AgentContainer();

    /// Add a single agent
    size_t add(agent_t* newAgent);
    /// Add all agents in the vector
    size_t add(std::vector<agent_t*> newAgents);

    /// Remove agent pointed to by targetPtr
    void remove(agent_t* targetPtr);
    /// Remove all agents of given typeId
    void remove(AgentTypeId typeId);
    /// Remove all agents for which the predicate p returns true
    template<class BinaryPred>void remove(BinaryPred p);
    /// Remove all agents from a target set
    void remove(std::set<agent_t*> targetSet);

    /// Replace one agent with another one in the same location, keeping bonds and neighbours
    /// old agent referenced by index, new one by a pointer
    /// new agent gets the old one's ID
    void replace(size_t iOldAgent, agent_t* newAgent);

    double getAverageStress() const;
    double getAverageStress(AgentTypeId typeId) const;
    double getAverageDeviation() const;
    double getMaxStress() const;

    double getMinCoord(size_t axis, AgentTypeId typeId = tAny) const;
    double getMaxCoord(size_t axis, AgentTypeId typeId = tAny) const;

    /// Count all agents
    size_t count() const;
    /// Count all agents of a given type
    size_t count(AgentTypeId typeId) const;
    size_t countSelected(AgentTypeId typeId) const;
    /// Count all agents for which the predicate p returns true
    template<class BinaryPred> size_t count(BinaryPred p) const;

    /// Apply functor f to all agents
    template<class BinaryFunction> void forAll(BinaryFunction f);
    /// Apply functor f to all agents of a given type
    template<class BinaryFunction> void forAll(AgentTypeId typeId, BinaryFunction f);
    /// Apply functor f to all agents for which the predicate p returns true
    template<class BinaryPred, class BinaryFunction> void forAll(BinaryPred p, BinaryFunction f);

    /// Synchronous agent rule execution.
    /// if agents are created or deleted by agent rules, they are added to or removed from the agent container
    /// once rules have been executed for all agents present at the time of calling this method.
    void execAgentRules();

    void writeFile(std::string const& fileName, AgentTypeId typeId=tAny) const;
    template<class BinaryPred> void writeFile(std::string const& fileName, BinaryPred p) const;

    /// Return a reference to the underlying data structure;
    /// this clearly breaks encapsulation, and should be used with great care.
    std::vector<agent_t*> & getAgentVector() { return agents; }
    std::vector<agent_t*> const& getAgentVector() const { return agents; }

private:
    // remove references to agent pointed to by targetPtr from all neighbourhood records.
    void removeReferences(agent_t *const targetPtr);
    // remove references to agent pointed to by targetPtr from its neighbours' records.
    void removeReferencesFromNeighbours(agent_t *const targetPtr);

    std::vector<agent_t*> agents;
};


} // end namespace absmc

#endif
