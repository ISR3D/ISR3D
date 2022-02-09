/*  This file is part of the Agent Based SMC Model for the COAST Project.
 *
 *  Copyright (C) 2008 Bernd Stahl
 *  Address: Battelle Batiment A, Route de Drize 7, 1227 Carouge, Switzerland
 *  E-mail: bernd.stahl@unige.ch
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef ABSMC_IMAGE_WRITER_HH
#define ABSMC_IMAGE_WRITER_HH

#include <cmath>
#include <limits>
#include <fstream>
#include <cstdlib>
#include <cstdio>

#include "graphics/imageWriter.h"

namespace absmc {

namespace graphics {

template<class Agent>
ImageWriter<Agent>::ImageWriter(double x0_, double y0_, double x1_, double y1_, double resolution_, std::string const& colorMapName)
    : x0(x0_), y0(y0_), x1(x1_), y1(y1_),
      resolution(resolution_),
      colorMap(0),
      imageTool(x0, y0, x1, y1, resolution),
      image(0)
{
    colorMap = generateColorMap(colorMapName) ;
    assert(colorMap);
    image = new PixMap<rgb>(imageTool.getNI(), imageTool.getNJ() );
    assert(image);
    startNew();
}

template<class Agent>
ImageWriter<Agent>::~ImageWriter() {
    delete colorMap;
    delete image;
}

template<class Agent>
void ImageWriter<Agent>::startNew()
{
  assert(image);
  image->setAll(rgb(0.0, 0.0, 0.0) );
}

template<class Agent>
void ImageWriter<Agent>::drawAgents(std::vector<agent_t*> const& agents,
                                    std::string const& quantity,
                                    double minVal, double maxVal)
{
    assert(image);

    if (fabs(maxVal-minVal)<1e-12) {
        computeMinMax(agents, quantity, minVal, maxVal);
    }

    const size_t numAgents = agents.size();
    for (size_t iAgent=0; iAgent<numAgents; iAgent++) {

        Agent const& agent = *(agents[iAgent]);

        // compute circle color
        double scaledValue = (agent.getScalarByName(quantity)-minVal) / (maxVal-minVal);
        scaledValue = std::max(0.0, scaledValue);
        scaledValue = std::min(1.0, scaledValue);

        rgb color = colorMap->get(scaledValue);

        // draw circle
        imageTool.drawCircle(*image, color, agent.getPos()[0], agent.getPos()[1], agent.getR(), true);
  }
}

template<class Agent>
void ImageWriter<Agent>::writeAsciiPpm(std::string const& fileName)
{
    assert(image);
    std::ofstream fout(fileName.c_str() );

    const size_t nI = imageTool.getNI();
    const size_t nJ = imageTool.getNJ();
    const int maxColor = 255;

   // write ppm header
    fout << "P3" << std::endl;
    fout << nI<< " " << nJ << std::endl;
    fout << maxColor << std::endl;

    // write image data, one rgb triplet per line
    for (size_t j=0; j<nJ; j++) {
        for (size_t i=0; i<nI; i++) {
            rgb const& color = image->getPixel(i, j);
            fout << static_cast<unsigned char>(color.r*maxColor) << " "
                 << static_cast<unsigned char>(color.g*maxColor) << " "
                 << static_cast<unsigned char>(color.b*maxColor) << " "
                 << std::endl;
        }
    }
}

template<class Agent>
void ImageWriter<Agent>::writePpm(std::string const& fileName)
{
    assert(image);

    const long unsigned int nI = (long unsigned int)imageTool.getNI();
    const long unsigned int nJ = (long unsigned int)imageTool.getNJ();
    const int maxColor = 255;

    // create memory-contiguous copy of image data
    unsigned char* data = new unsigned char[nI*nJ*3];
    for (size_t i=0; i<nI; i++) {
        for (size_t j=0; j<nJ; j++) {
            rgb const& color = image->getPixel(i, j);
            data[(j*nI+i)*3+0] = static_cast<unsigned char>(color.r*maxColor);
            data[(j*nI+i)*3+1] = static_cast<unsigned char>(color.g*maxColor);
            data[(j*nI+i)*3+2] = static_cast<unsigned char>(color.b*maxColor);
        }
    }


    // write binary ppm file
    FILE* fp = fopen(fileName.c_str(), "wb");
    fprintf(fp, "P6\n%lu %lu\n%d\n", nI, nJ, maxColor);
    fwrite(data, sizeof(unsigned char), nI*nJ*3, fp);
    fclose(fp);

    delete data;
}


template<class Agent>
void ImageWriter<Agent>::writePng(std::string const& fileName)
{
    writePpm(fileName+".ppm");

    std::string convCommand =
            std::string("convert ") + fileName + ".ppm " + fileName + ".png";
    std::string rmCommand =
            std::string("/bin/rm ") + fileName + ".ppm";
    system(convCommand.c_str());
    system(rmCommand.c_str());
}

//////////////////// helper functions ///////////////////////////////////////

template<class Agent>
void computeMinMax(std::vector<Agent*> const& agents,
                   std::string const& quantity,
                   double& minVal, double& maxVal)
{
    minVal =  std::numeric_limits<double>::max();
    maxVal = -std::numeric_limits<double>::max();

    const size_t numAgents = agents.size();
    for(size_t iAgent=0; iAgent<numAgents; iAgent++) {
        Agent const& agent = *(agents[iAgent]);
        minVal = std::min(minVal, agent.getScalarByName(quantity) );
        maxVal = std::max(maxVal, agent.getScalarByName(quantity) );
    }
    std::cout << "ImageWriter auto-scale: minValue=" << minVal << "  maxValue=" << maxVal << std::endl;
}

} // end namespace graphics

} // end namespace absmc

#endif
