#ifndef ABSMC_AGENT_CONTAINER_HH
#define ABSMC_AGENT_CONTAINER_HH

#include <vector>
#include <cassert>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <limits>

#include "core/agentContainer.h"


namespace absmc {

template<size_t nDim>
AgentContainer<nDim>::~AgentContainer()
{
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        delete agents[iAgent];
    }
}

template<size_t nDim>
size_t AgentContainer<nDim>::add(agent_t* newAgent)
{
    agents.push_back(newAgent);
    return agents.size()-1;
}


template<size_t nDim>
size_t AgentContainer<nDim>::add(std::vector<agent_t*> newAgents)
{
    std::copy(newAgents.begin(), newAgents.end(), std::back_inserter(agents) );
    return agents.size()-1;
}

template<size_t nDim>
void AgentContainer<nDim>::removeReferences(agent_t *const targetPtr)
{
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        auto& neighbours = agents[iAgent]->getNeighbours();
        neighbours.erase(std::remove(neighbours.begin(), neighbours.end(), targetPtr), neighbours.end() );

        auto& bonds = agents[iAgent]->getBonds();
        bonds.erase(std::remove(bonds.begin(), bonds.end(), targetPtr), bonds.end() );
    }
}

// template<size_t nDim>
// void AgentContainer<nDim>::removeReferences(std::set<agent_t*> targetSet)
// {
//     const size_t nAgents = agents.size();
//     for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
//         typename AgentBase<nDim>::nb_vec_t & neighbours = agents[iAgent]->getNeighbours();
//         for (auto& agentPtr : neighbours) {
//             if (targetSet.find(agentPtr) != targetSet.end())
//         }
//         neighbours.erase(std::remove(neighbours.begin(), neighbours.end(), targetPtr), neighbours.end() );
//     }
//
//     for (auto agentPtr : agents) {
//         if (targetSet.find(agentPtr) != targetSet.end()) {
//             delete agentPtr;
//             agentPtr = nullptr;
//         }
//     }
//
//     agents.erase(std::remove(agents.begin(), agents.end(), nullptr), agents.end() );
// }

// we assume here that the neighbourhoods and bonds are symmetric and valid
template<size_t nDim>
void AgentContainer<nDim>::removeReferencesFromNeighbours(agent_t *const targetPtr)
{
    auto neighbourList = targetPtr->getNeighbours();
    for (auto& nbr : neighbourList) {
        auto& neighbours = nbr->getNeighbours();
        neighbours.erase(std::remove(neighbours.begin(), neighbours.end(), targetPtr), neighbours.end() );
    }
    auto bondList = targetPtr->getBonds();
    for (auto& bnd : bondList) {
        bnd.other->removeBond(targetPtr);
    }
}

template<size_t nDim>
void AgentContainer<nDim>::remove(agent_t* targetPtr)
{
    // remove references to agent pointed to by targetPtr from all neighbourhood records.
    removeReferences(targetPtr);
    // delete target agent and remove pointer from agent vector
    delete targetPtr;
    agents.erase(std::remove(agents.begin(), agents.end(), targetPtr), agents.end() );
}

template<size_t nDim>
void AgentContainer<nDim>::remove(std::set<agent_t*> targetSet)
{

    /// openMP version has to be checked before being used; so far, single-threaded version works sufficiently fast
    /// remove references to agent pointed to by targetPtr from all neighbourhood records.
    // #pragma omp parallel for schedule(dynamic, 10000)
    // for (auto it = agents.begin(); it < agents.end(); ++it) {
    //     if (targetSet.find(*it) != targetSet.end()) {
    //         removeReferencesFromNeighbours(*it);
    //     }
    // }
    for (auto& agentPtr : agents) {
        if (targetSet.find(agentPtr) != targetSet.end()) {
            removeReferencesFromNeighbours(agentPtr);
        }
    }

    // #pragma omp parallel for schedule(dynamic, 10000)
    // for (auto it = agents.begin(); it < agents.end(); ++it) {
    //     if (targetSet.find(*it) != targetSet.end()) {
    //         delete *it;
    //         *it = nullptr;
    //     }
    // }
    for (auto& agentPtr : agents) {
        if (targetSet.find(agentPtr) != targetSet.end()) {
            delete agentPtr;
            agentPtr = nullptr;
        }
    }

    agents.erase(std::remove(agents.begin(), agents.end(), nullptr), agents.end() );
}

template<size_t nDim>
void AgentContainer<nDim>::remove(AgentTypeId typeId)
{
    const size_t nAgents = agents.size();

    // first loop: remove reference,
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (agents[iAgent]->getTypeId()==typeId) {
            removeReferencesFromNeighbours(agents[iAgent]);
        }
    }
    // second loop: delete target agents, set corresponding pointers to nil
    agent_t *const nil = 0;
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (agents[iAgent]->getTypeId()==typeId) {
            delete agents[iAgent];
            agents[iAgent] = nil;
        }
    }
    // third loop: remove all nil pointers
    agents.erase(std::remove(agents.begin(), agents.end(), nil), agents.end() );
}

template<size_t nDim>
template<class BinaryPred>
void AgentContainer<nDim>::remove(BinaryPred p)
{
    const size_t nAgents = agents.size();
    std::set<agent_t*> setToRemove;
    // first loop: remove reference,
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (p(iAgent, agents[iAgent]) ) {
            setToRemove.insert(agents[iAgent]);
        }
    }
    remove(setToRemove);
}


/// Replace one agent with another one in the same location, keeping bonds and neighbours
template<size_t nDim>
void AgentContainer<nDim>::replace(size_t iOldAgent, agent_t* newAgent) {
    agent_t* oldAgent = agents[iOldAgent];

    newAgent->id            = oldAgent->id;
    newAgent->pos           = oldAgent->pos;
    newAgent->pos0          = oldAgent->pos0;
    newAgent->neighbours    = oldAgent->neighbours;
    newAgent->bonds         = oldAgent->bonds;
    if (oldAgent->surfaceNormal != nullptr)
        newAgent->surfaceNormal = new point_t(*(oldAgent->surfaceNormal));
    newAgent->isMobile      = oldAgent->isMobile;
    newAgent->isAffected    = oldAgent->isAffected;

    auto neighbourList = newAgent->getNeighbours();
    for (auto& nbr : neighbourList) {
        auto& neighbours = nbr->getNeighbours();
        for (auto& backLink : neighbours) {
            if(backLink == oldAgent) {
                backLink = newAgent;
            }
        }
    }
    auto bondList = newAgent->getBonds();
    for (auto& bnd : bondList) {
        auto& bonds = bnd.other->getBonds();
        for (auto& backLink : bonds) {
            if(backLink.other == oldAgent) {
                backLink.other = newAgent;
            }
        }
    }

    agents[iOldAgent] = newAgent;
    delete oldAgent;
}


/// general method "get average quantity" should be implemented if these methods are to be used
template<size_t nDim>
double AgentContainer<nDim>::getAverageStress(AgentTypeId typeId) const
{
    double strs =  0.0;
    int count = 0;
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (agents[iAgent]->getTypeId() == typeId) {strs += norm<3>(agents[iAgent]->getStress()); count++;}
    }
    strs = strs/count;
    return strs;
}

template<size_t nDim>
double AgentContainer<nDim>::getAverageDeviation() const
{
    double strs =  0.0;
    int count = 0;
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        strs += fabs(norm<3>(agents[iAgent]->getPos())-norm<3>(agents[iAgent]->getPos0())); count++;
    }
    strs = strs/count;
    return strs;
}

template<size_t nDim>
double AgentContainer<nDim>::getAverageStress() const
{
    double strs =  0.0;
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        strs += norm<3>(agents[iAgent]->getStress());
    }
    strs = strs/agents.size();
    return strs;
}

template<size_t nDim>
double AgentContainer<nDim>::getMaxStress() const
{
    double strs =  0.0;
    double strsmax =  0.0;
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        strs = norm<3>(agents[iAgent]->getStress());
        if (strs > strsmax) {strsmax = strs;}
    }

    return strsmax;
}

/// gets the max coordinate along the i-th axis
template<size_t nDim>
double AgentContainer<nDim>::getMaxCoord(size_t axis, AgentTypeId typeId) const
{
    assert(axis <= 2);
    double coordmax = std::numeric_limits<double>::lowest();
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (typeId == tAny || typeId == agents[iAgent]->getTypeId()){
            double cur = agents[iAgent]->getPos()[axis];
            if (cur > coordmax) {coordmax = cur;}
        }
    }

    return coordmax;
}

/// gets the min coordinate along the i-th axis
template<size_t nDim>
double AgentContainer<nDim>::getMinCoord(size_t axis, AgentTypeId typeId) const
{
    assert(axis <= 2);
    double coordmin = std::numeric_limits<double>::max();
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (typeId == tAny || typeId == agents[iAgent]->getTypeId()){
            double cur = agents[iAgent]->getPos()[axis];
            if (cur < coordmin) {coordmin = cur;}
        }
    }

    return coordmin;
}

template<size_t nDim>
inline size_t AgentContainer<nDim>::count() const
{
    return agents.size();
}

template<size_t nDim>
size_t AgentContainer<nDim>::count(AgentTypeId typeId) const
{
    size_t nMatching = 0;
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (agents[iAgent]->getTypeId() == typeId) nMatching++;
    }
    return nMatching;
}

template<size_t nDim>
template<class BinaryPred>
size_t AgentContainer<nDim>::count(BinaryPred p) const
{
    size_t nMatching = 0;
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (p(iAgent, agents[iAgent]) ) nMatching++;
    }
    return nMatching;
}

template<size_t nDim>
template<class BinaryFunction>
void AgentContainer<nDim>::forAll(BinaryFunction f)
{
    const size_t nAgents = agents.size();
    #pragma omp parallel for schedule(dynamic, 10000)
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        f(iAgent, agents[iAgent]);
    }
}

template<size_t nDim>
template<class BinaryFunction>
void AgentContainer<nDim>::forAll(AgentTypeId typeId, BinaryFunction f)
{
    const size_t nAgents = agents.size();
    #pragma omp parallel for schedule(dynamic, 10000)
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (agents[iAgent]->getTypeId() == typeId) f(iAgent, agents[iAgent]);
    }
}


template<size_t nDim>
template<class BinaryPred, class BinaryFunction>
void AgentContainer<nDim>::forAll(BinaryPred p, BinaryFunction f)
{
    const size_t nAgents = agents.size();
    #pragma omp parallel for schedule(dynamic, 10000)
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (p(iAgent, agents[iAgent]) ) f(iAgent, agents[iAgent]);
    }
}

template<size_t nDim>
void AgentContainer<nDim>::execAgentRules()
{
    /// If new agents are created by agent rules, rules are not executed on these agents in this loop,
    /// since the upper loop bound is determined before; this is exactly the desired behaviour.
    /// If agents are deleted by agent rules, they are removed from the agent container only after completing the loop.
    std::vector<agent_t*> agentsCreated, agentsDeleted;
    std::vector<agent_t*> agentsToRemove;

    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        agentsCreated.clear();
        agentsDeleted.clear();
        agents[iAgent]->execAgentRule(iAgent, agentsCreated, agentsDeleted);
        // agents created can simply be added, see above.
        std::copy(agentsCreated.begin(), agentsCreated.end(), std::back_inserter(agents) );
        // agents deleted must be removed afer loop completion.
        std::copy(agentsDeleted.begin(), agentsDeleted.end(), std::back_inserter(agentsToRemove) );
    }

    /*const size_t nAgentsToRemove = agentsToRemove.size();
    for (size_t iAgent=0; iAgent<nAgentsToRemove; iAgent++) {
        remove(agentsToRemove[iAgent]);
    }*/

    std::set<agent_t*> setToRemove(agentsToRemove.begin(), agentsToRemove.end());
    remove(setToRemove);
}

template<size_t nDim>
void AgentContainer<nDim>::writeFile(std::string const& fileName, AgentTypeId typeId) const
{
    std::ofstream ofstr (fileName.c_str());
    if (!ofstr) {
        std::cout << "Cannot open file " << fileName << " for writing." << std::endl;
        return;
    }
    ofstr << std::setprecision(8);

    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (typeId==tAny || agents[iAgent]->getTypeId()==typeId) {
            agents[iAgent]->print(ofstr);
            ofstr << std::endl;
        }
    }
}

template<size_t nDim>
template<class BinaryPred>
void AgentContainer<nDim>::writeFile(std::string const& fileName, BinaryPred p) const
{
    std::ofstream ofstr (fileName.c_str());
    if (!ofstr) {
        std::cout << "Cannot open file " << fileName << " for writing." << std::endl;
        return;
    }
    ofstr << std::setprecision(8);

    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        if (p(iAgent, agents[iAgent]) ) {
            agents[iAgent]->print(ofstr);
            ofstr << std::endl;
        }
    }
}

} // end namespace absmc

#endif
