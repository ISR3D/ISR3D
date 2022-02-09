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

#ifndef ABSMC_IMAGE_WRITER_H
#define ABSMC_IMAGE_WRITER_H

#include <string>
#include <vector>

#include "graphics/colormaps.h"
#include "graphics/imageTool.h"
#include "graphics/pixMap.h"

namespace absmc {

namespace graphics {

template<class Agent>
class ImageWriter {
public:
    typedef Agent agent_t;

    ImageWriter(double x0_, double y0_, double x1_, double y1_, double resolution_, std::string const& colorMapName="default");

    ~ImageWriter();

    void startNew();

    /// draw agents represented by circles, colored according to the given quantity and color map
    void drawAgents(std::vector<agent_t*> const& agents,
                    std::string const& quantity,
                    double minVal, double maxVal);

    /// write ascii ppm file; too slow in practice
    void writeAsciiPpm(std::string const& fileName);
    /// write binary ppm file; platform-indepence is uncertain
    void writePpm(std::string const& fileName);
    /// write binary ppm file make png file using the 'convert' program
    void writePng(std::string const& fileName);

private:
    const double x0, y0, x1, y1;            // physical coordinates       [simulation_units]
    const double resolution;                 // in pixels/physical unit    [1/simulation_units]
    ColorMap* colorMap;
    ImageTool<rgb> imageTool;
    PixMap<rgb>* image;
};

template<class Agent>
void computeMinMax(std::vector<Agent*> const& agents,
                   std::string const& quantity,
                   double& minVal, double& maxVal);

} // end namespace graphics

} // end namespace absmc

#endif
