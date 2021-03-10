#ifndef ABSMC_AGENT_FILE_READER_HH
#define ABSMC_AGENT_FILE_READER_HH

#include <fstream>
#include <cassert>
#include <string>

#include "core/agentFileReader.h"
//#include "util/miscUtil.h"

namespace absmc {


template<size_t nDim>
size_t AgentFileReader<nDim>::readFile(std::string const& fileName,
                                       AgentFactory<nDim> const& factory,
                                       std::vector<agent_t*> & agents,
                                       AgentTypeId typeId)
{
    std::ifstream ifstr(fileName.c_str());
    if (!ifstr) {
        std::cout << "Cannot open file " << fileName << " for reading." << std::endl;
        return 0;
    }
    size_t numCreated = 0;

    AgentTypeId tId;
    point_t pos, pos0;
    double r;
    double drugConc, wssOsi, wssMax, neoDist = 0, nitricConc = 0;
    int age, clock, lengthG1;
    bool baFlag, ciFlag, selectedForNOFlag;
    while (ifstr) {
        ifstr >> tId;
        for (size_t iDim=0; iDim<nDim; iDim++) ifstr >> pos[iDim];
        for (size_t iDim=0; iDim<nDim; iDim++) ifstr >> pos0[iDim];
        ifstr >> r;
        ifstr >> drugConc;
        ifstr >> wssOsi;
        ifstr >> wssMax;
        if (tId != tObstacle3D) {
            ifstr >> neoDist; //TODO: get back from recent stage3 file
            if (tId == tSMC3D) {
                ifstr >> age;
                ifstr >> clock;
                ifstr >> lengthG1;
                ifstr >> baFlag;
                ifstr >> ciFlag;
                //ifstr >> pmFlag;
                ifstr >> selectedForNOFlag; // TODO: set and match with nitricConc
            }
            if (tId == tEndo3D) {
                ifstr >> age;
                ifstr >> clock;
                ifstr >> lengthG1;
                ifstr >> baFlag;
                ifstr >> ciFlag;
            }
        }
        if (typeId==tAny || tId==typeId) {
            agent_t* newAgent = (ifstr) ? factory.create(tId, pos, pos0, r, drugConc, wssOsi, wssMax, neoDist, nitricConc,
                                                         age, clock, lengthG1,
                                                         baFlag, ciFlag, false) : 0;
            if (newAgent) {
                agents.push_back(newAgent);
                numCreated++;
            }
        }
    }
    return numCreated;
}




template<size_t nDim>
size_t AgentFileReader<nDim>::readFile(std::string const& fileName,
                                       AgentFactory<nDim> const& factory,
                                       AgentContainer<nDim> & agentContainer,
                                       AgentTypeId typeId)
{
    std::vector<agent_t*> agents;
    readFile(fileName, factory, agents, typeId);
    agentContainer.add(agents);
    return agents.size();
}




template<size_t nDim>
size_t AgentFileReader<nDim>::readCsvFile(std::string const& fileName,
                                       AgentFactory<nDim> const& factory,
                                       std::vector<agent_t*> & agents)
{
    std::ifstream ifstr(fileName.c_str());
    if (!ifstr) {
        std::cout << "Cannot open file " << fileName << " for reading." << std::endl;
        return 0;
    }
    size_t numCreated = 0;

    AgentTypeId tId;
    point_t pos, pos0;
    double r;
    // double drugConc, wssOsi, wssMax, neoDist = 0, nitricConc = 0;
    // int age, clock, lengthG1;
    // bool baFlag, ciFlag, selectedForNOFlag;
    while (ifstr) {
        ifstr >> tId;

        for (size_t iDim=0; iDim<nDim; iDim++) ifstr >> pos[iDim];
        pos0 = pos;
        ifstr >> r;

        // ifstr >> drugConc;
        // ifstr >> wssOsi;
        // ifstr >> wssMax;
        // if (tId != tObstacle3D) {
        //     ifstr >> neoDist; //TODO: get back from recent stage3 file
        //     if (tId == tSMC3D) {
        //         ifstr >> age;
        //         ifstr >> clock;
        //         ifstr >> lengthG1;
        //         ifstr >> baFlag;
        //         ifstr >> ciFlag;
        //         //ifstr >> pmFlag;
        //         ifstr >> selectedForNOFlag; // TODO: set and match with nitricConc
        //     }
        //     if (tId == tEndo3D) {
        //         ifstr >> age;
        //         ifstr >> clock;
        //         ifstr >> lengthG1;
        //         ifstr >> baFlag;
        //         ifstr >> ciFlag;
        //     }
        // }
        // agent_t* newAgent = (ifstr) ? factory.create(tId, pos, pos0, r, drugConc, wssOsi, wssMax, neoDist, nitricConc,
        //                                              age, clock, lengthG1,
        //                                              baFlag, ciFlag, false) : 0;

        agent_t* newAgent = (ifstr) ? factory.create(tId, pos, pos0, r) : 0;

        if (newAgent) {
            agents.push_back(newAgent);
            numCreated++;
        }
    }
    return numCreated;
}




template<size_t nDim>
size_t AgentFileReader<nDim>::readCsvFile(std::string const& fileName,
                                       AgentFactory<nDim> const& factory,
                                       AgentContainer<nDim> & agentContainer)
{
    std::vector<agent_t*> agents;
    readCsvFile(fileName, factory, agents);
    agentContainer.add(agents);
    return agents.size();
}

} // end namespace absmc

#endif
