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
 
#ifndef ABSMC_IMAGE_TOOL_H
#define ABSMC_IMAGE_TOOL_H

#include "graphics/rgb.h"
#include "graphics/pixMap.h"

namespace absmc {

namespace graphics {

template<typename T>
class ImageTool {
public:
    ImageTool(double x0_, double y0_, double x1_, double y1_, double res_);
    
    size_t getNI() const  { return nI; }
    size_t getNJ() const  { return nJ; }
    
    size_t getI(double x) const;
    size_t getJ(double y) const;
    
    void drawLine(PixMap<T> & pixMap, T const& pixelValue,
                  double xa, double ya, double xb, double yb) const;
    
    void drawCircle(PixMap<T> & pixMap, T const& pixelValue,
                   double x, double y, double r, bool fill=true) const;
private:
    const double x0, y0, x1, y1, res;
    const size_t nI, nJ;
};

} // end namespace graphics

} // end namespace absmc

#endif
