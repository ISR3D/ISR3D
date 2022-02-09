#ifndef ABSMC_ANALYSIS_UTILS_3D_H
#define ABSMC_ANALYSIS_UTILS_3D_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>

#include "core/geometry.h"
#include "core/euclideanMetrics.h"
#include "core/agentBase.h"
#include "core/agentContainer.h"

namespace absmc {

/// For every agent in container, print complete neighbourhood record (one agent per line).
void printNB(AgentContainer<3> & agentContainer, std::string const& fileName)
{
    const size_t nDim = 3;
    typedef AgentBase<nDim> agent_t;
    //typedef Point<nDim, double> point_t;

    std::ofstream ofstr(fileName.c_str() );
    if (!ofstr) {
        std::cout << "Can't open file " << fileName << " for writing." << std::endl;
        return;
    }

    std::vector<agent_t*> const& agents = agentContainer.getAgentVector();
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        agent_t::nb_vec_t const& neighbours = agents[iAgent]->getNeighbours();

        for (size_t iNb=0; iNb<neighbours.size(); iNb++) {
            ofstr << neighbours[iNb] << "  ";
        }
        ofstr << std::endl;
    }
}

/// For every agent in container, print distances to every neighbour: dX, dY, dZ, r, pair of agents per line.
void printNBDistances(AgentContainer<3> & agentContainer, std::string const& fileName)
{
    const size_t nDim = 3;
    typedef AgentBase<nDim> agent_t;
    typedef Point<nDim, double> point_t;
    typedef EuclideanMetrics<3, double> metrics_t;

    std::ofstream ofstr(fileName.c_str() );
    if (!ofstr) {
        std::cout << "Can't open file " << fileName << " for writing." << std::endl;
        return;
    }

    metrics_t metrics;

    std::vector<agent_t*> const& agents = agentContainer.getAgentVector();
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {

        point_t const& pos1 = agents[iAgent]->getPos();

        agent_t::nb_vec_t const& neighbours = agents[iAgent]->getNeighbours();
        for (size_t iNb=0; iNb<neighbours.size(); iNb++) {

            point_t const& pos2 = neighbours[iNb]->getPos();

            const double dX = pos2[0]-pos1[0];
            const double dY = pos2[1]-pos1[1];
            const double dZ = pos2[2]-pos1[2];
            const double r = metrics.dist(pos1, pos2);
            ofstr << dX << " " << dY << " " << dZ << " " << r << std::endl;
        }
    }
}

/// For every agent in container, print number of neighbours (one agent per line).
void printNBCount(AgentContainer<3> & agentContainer, std::string const& fileName)
{
    const size_t nDim = 3;
    typedef AgentBase<nDim> agent_t;
    typedef Point<nDim, double> point_t;

    std::ofstream ofstr(fileName.c_str() );
    if (!ofstr) {
        std::cout << "Can't open file " << fileName << " for writing." << std::endl;
        return;
    }

    std::vector<agent_t*> const& agents = agentContainer.getAgentVector();
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        point_t const& pos = agents[iAgent]->getPos();
        agent_t::nb_vec_t const& neighbours = agents[iAgent]->getNeighbours();
        ofstr << pos << " " << neighbours.size() << " " << std::endl;
    }
}


/// Compute mean and standard deviation of the distribution of distance-to-neighbours.
template<class Agent>
void computeNbDistStat(AgentContainer<3> & agentContainer, double* mean, double* stdDev)
{
    const size_t nDim = 3;
    typedef AgentBase<nDim> agent_t;
    typedef Point<nDim, double> point_t;
    typedef EuclideanMetrics<Agent::nDim, double> metrics_t;

    metrics_t metrics;

    std::vector<double> valX;

    std::vector<agent_t*> const& agents = agentContainer.getAgentVector();
    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {

        point_t const& pos1 = agents[iAgent]->getPos();
        agent_t::nb_vec_t const& neighbours = agents[iAgent]->getNeighbours();

        for (size_t iNb=0; iNb<neighbours.size(); iNb++) {

            point_t const& pos2 = neighbours[iNb]->getPos();
            const double r = metrics.dist(pos1, pos2);
            valX.push_back(r);
        }
    }

    const size_t N = valX.size();
    (*mean) = std::accumulate(valX.begin(), valX.end(), 0.0) / N;

    double sumDiffSqr = 0.0;
    for (size_t i=0; i<N; i++) {
        sumDiffSqr += (valX[i]-(*mean))*(valX[i]-(*mean));
    }

    (*stdDev) = sumDiffSqr/N;
}

/// Compute l2 norm of difference between initial and current positions.
template<class Agent>
double l2Difference(std::vector<Agent*> const& agents)
{
    typedef Point<Agent::nDim, double> point_t;

    double sumSquares = 0.0;

    const size_t nAgents = agents.size();
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {

        const point_t pos = agents[iAgent]->getPos();
        const point_t pos0 = agents[iAgent]->getPos0();

        point_t deltaPos;
        for (size_t iDim=0; iDim<Agent::nDim; iDim++) {
            deltaPos[iDim] = pos[iDim] - pos0[iDim];
        }
        sumSquares += normSqr<Agent::nDim>(deltaPos);
    }
    return sqrt(sumSquares / nAgents);
}

} // end namespace absmc

#endif
