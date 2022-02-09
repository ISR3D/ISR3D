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
 
#ifndef ABSMC_IMAGE_TOOL_HH
#define ABSMC_IMAGE_TOOL_HH

#include <cmath>
#include <algorithm>

#include "graphics/imageTool.h"
#include "util/mathUtil.h"

namespace absmc {

namespace graphics {

template<typename T>
ImageTool<T>::ImageTool(double x0_, double y0_, double x1_, double y1_, double res_)
    : x0(x0_), y0(y0_), x1(x1_), y1(y1_), res(res_),
      nI(static_cast<size_t>(ceil( fabs(x1-x0) * res)) ),
      nJ(static_cast<size_t>(ceil( fabs(y1-y0) * res)) )
{ }
      
template<typename T>
inline size_t ImageTool<T>::getI(double x) const
{
    size_t i = static_cast<size_t>((x-x0)*res);
    i = std::max((size_t)0, i);
    i = std::min(nI-1, i);
    return i;
}

template<typename T>
inline size_t ImageTool<T>::getJ(double y) const
{
    // flip j coordinate such that y increases from bottom to top
    size_t j = (nJ-1) - static_cast<size_t>((y-y0)*res);
    j = std::max((size_t)0, j);
    j = std::min(nJ-1, j);
    return j;
}

template<typename T>
void ImageTool<T>::drawLine(PixMap<T> & pixMap, T const& pixelValue,
                         double xa, double ya, double xb, double yb) const
{
    const size_t numSteps = (xb!=xa) ?
            static_cast<size_t>(fabs(xb-xa)*res) : static_cast<size_t>(fabs(yb-ya)*res);

    const double dX = (xb-xa) / (numSteps-1);
    const double dY = (yb-ya) / (numSteps-1);

    for (size_t iStep=0; iStep<numSteps; iStep++) {
        size_t i = getI(xa + iStep*dX);
        size_t j = getJ(ya + iStep*dY);
        pixMap.setPixel(i, j, pixelValue);
    }
}

template<typename T>
void ImageTool<T>::drawCircle(PixMap<T> & pixMap, T const& pixelValue,
                           double x, double y, double r, bool fill) const
{
    // calculate bounding box
    size_t iLo = getI(x-r);
    size_t iHi = getI(x+r);
    // flip j coordinate such that y increases from bottom to top
    size_t jLo = getJ(y+r);
    size_t jHi = getJ(y-r);
    
    const double rSqr = r*r;
    // loop over pixels in bounding box only
    for (size_t i=iLo; i<=iHi; i++ ) {
        for (size_t j=jLo; j<=jHi; j++ ) {
            const double xDiff = ((double)i/res + x0) - x;
            const double yDiff = ((double)(nJ-1-j)/res + y0) - y;
            const double dSqr = xDiff*xDiff + yDiff*yDiff;

            if (dSqr <= rSqr) {
                // draw filled circle
                if (fill) {
                    pixMap.setPixel(i, j, pixelValue);
                }
                // draw circle outline
                else if (dSqr>rSqr-4*r/res+4/res/res) {
                    pixMap.setPixel(i, j, pixelValue);
                }
            }
        }
    }
}


} // end namespace graphics

} // end namespace absmc

#endif
