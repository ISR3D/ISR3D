#ifndef ABSMC_VTK_WRITER_HH
#define ABSMC_VTK_WRITER_HH

#include <string>
#include <fstream>
#include <sstream>
#include <limits>

#include "graphics/vtkWriter.h"
#include <util/fileUtil.h>

namespace absmc {

namespace graphics {

template<class Agent>
void VtkWriter<Agent>::writeVtkPolyData(std::vector<agent_t*> const& agents,
                                        std::string const& fileName,
                                        std::string const& scalarQuantityString,
                                        std::string const& vectorQuantityString)
{
    const size_t nAgents = agents.size();

    // parse quantity strings
    std::vector<std::string> scalarQuantities;
    util::tokenizeString(scalarQuantityString, scalarQuantities, '/');
    if (!scalarQuantities.empty() && scalarQuantities[0]=="none") { scalarQuantities.clear(); }
    std::vector<std::string> vectorQuantities;
    util::tokenizeString(vectorQuantityString, vectorQuantities, '/');
    if (!vectorQuantities.empty() && vectorQuantities[0]=="none") { vectorQuantities.clear(); }

    // open fstream for writing
    std::ofstream ofstr(fileName.c_str());
    if (!ofstr) {
        std::cout << "Cannot open file " << fileName << " for writing." << std::endl;
        return;
    }

    // write vtp header
    ofstr << "<?xml version=\"1.0\"?>" << std::endl;
    ofstr <<  "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    ofstr << "<PolyData>" << std::endl;
    ofstr << "<Piece NumberOfPoints=\"" << nAgents<< "\" ";
    ofstr << "NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">" << std::endl;

    // write agent center coordinates
    ofstr << "<Points>" << std::endl;
    ofstr << "<DataArray type=\"Float64\" Name=\"position\" NumberOfComponents=\"3\">" << std::endl;
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        Agent * agent = agents[iAgent];
        const point_t pos = agent->getPos();
        ofstr << (float) pos[0] << " " << (float) pos[1] << " ";
        if (Agent::nDim>2) ofstr << (float) pos[2] << std::endl;
        else ofstr << 0 << std::endl;
    }
    ofstr << "</DataArray>" << std::endl;
    ofstr << "</Points>" << std::endl;

    // write agent radii
    ofstr << "<PointData Scalars=\"radius\">" << std::endl;
    ofstr << "<DataArray type=\"Float64\" Name=\"radius\">" << std::endl;
    for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
        Agent * agent = agents[iAgent];
        ofstr << (float) agent->getR() << " ";
    }
    ofstr << std::endl;
    ofstr << "</DataArray>" << std::endl;

    // write additional per-agent data as requested by user
    // scalar quantities
    const size_t nSQ = scalarQuantities.size();
    for (size_t iSQ=0; iSQ<nSQ; iSQ++) {
        const std::string quantity = scalarQuantities[iSQ];
        ofstr << "<DataArray type=\"Float64\" Name=\"" <<  quantity << "\">" << std::endl;
        for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
            Agent * agent = agents[iAgent];
            ofstr << (float) convertNaNRep(agent->getScalarByName(quantity)) << " ";
        }
        ofstr << std::endl;
        ofstr << "</DataArray>" << std::endl;
    }
    // vector quantities
    const size_t nVQ = vectorQuantities.size();
    for (size_t iVQ=0; iVQ<nVQ; iVQ++) {
        const std::string quantity = vectorQuantities[iVQ];
        ofstr << "<DataArray type=\"Float64\" Name=\"" <<  quantity << "\" NumberOfComponents=\"3\">" << std::endl;
        for (size_t iAgent=0; iAgent<nAgents; iAgent++) {
            Agent * agent = agents[iAgent];
            ofstr << convertNaNRep(agent->getVectorByName(quantity)) << " ";
        }
        ofstr << std::endl;
        ofstr << "</DataArray>" << std::endl;
    }

    ofstr << "</PointData>" << std::endl;

    // write closing xml tags
    ofstr << "</Piece>" << std::endl;
    ofstr << "</PolyData>" << std::endl;
    ofstr << "</VTKFile>" << std::endl;
}

template<class Agent>
inline double VtkWriter<Agent>::convertNaNRep(double scalar, double NaNRep)
{
    return (std::isnan(scalar)) ? NaNRep : scalar;
}

template<class Agent>
inline typename VtkWriter<Agent>::point_t VtkWriter<Agent>::convertNaNRep(point_t const& vector, double NaNRep)
{
    point_t v(vector);
    for (size_t iDim=0; iDim<Agent::nDim; iDim++) {
        if (std::isnan(v[iDim])) v[iDim] = NaNRep;
    }
    return v;
}


} // end namespace graphics

} // end namespace absmc

#endif
