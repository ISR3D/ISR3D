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
 
#include <cmath>
// #include <iostream>
 
#include "model2D/geometryGenerator2D.h"
#include "util/mathUtil.h"
#include "util/random.h"

namespace absmc {
    
size_t GeometryGenerator2D::generateHexPack(AgentGeoVector &  coordinates,
                                          double x0, double y0,                    
                                          size_t nX, size_t nY, double R,
                                          size_t nL)
{
    const double s  = sqrt(3.0)*R;
    
    double x, y;
    size_t numAgents = 0;
    for (size_t j=1; j<=nY; j++) {
        
        if (j % 2 == 0) {
            // Rod's code: i=1..nX+1
            // my suggestion: i=2..nX for cell-like boundaries, i=1..nX for true periodic boundaries
            for (size_t i=2; i<=nX; i++) {
                x = x0 + (i-1) * s;
                y = y0 + (3*j-1) * 0.5 * R;
                coordinates.push_back(AgentGeometry2D(x, y, x, y, R) );
                numAgents++;
            }
        } 
        else {
            for (size_t i=1; i<=nX; i++) {
                x = x0 + (i-0.5) * s;
                y = y0 + (3*j-1) * 0.5 * R;
                coordinates.push_back(AgentGeometry2D(x, y, x, y, R) );
                numAgents++;
            }
        }
    }
    return numAgents;
}
    
    
size_t GeometryGenerator2D::generateLiningLayer(AgentGeoVector &  liningCoordinates,
                                              AgentGeoVector &  cellCoordinates,
                                              size_t nX, size_t nL)
{
    if (nL<2 || nL>5) return 0;
    
    // coefficients as computed by Rod; cf theory document
    const double rLCoeff[6] = { 0.0, 0.0, 0.533, 0.281, 0.192, 0.146 };
    const double rACoeff[6] = { 0.0, 0.0, 1.3416, 1.1513, 1.0914, 1.0632 };
    
    size_t numLining = 0;
    for (size_t iAgent=0; iAgent<nX; iAgent++) {
        const double x = cellCoordinates[iAgent].pos[0]; 
        const double y = cellCoordinates[iAgent].pos[1];
        const double R = cellCoordinates[iAgent].r;
        
        const double rL = (rLCoeff[nL])*R;
        const double rA = (rACoeff[nL])*R;
        const double s  = sqrt(3.0)*R;
        const double theta0 = 0.5*math::pi - acos(0.5*s/(rA));
        const double dTheta = theta0 / (nL-1.0);
        
        const size_t numL = (iAgent==nX-1) ? 2*nL-1 : 2*nL-2;
        for (size_t iL=0; iL<numL; iL++) {
            double thetaN = -theta0 + iL*dTheta;
            double xL = x + rA*sin(thetaN);
            double yL = y - rA*cos(thetaN);
            liningCoordinates.push_back(AgentGeometry2D(xL, yL, xL, yL, rL) );
            numLining++;
        }        
    }
    return numLining;       
}
    
    
size_t GeometryGenerator2D::generatePerturbedHexPack(AgentGeoVector &  coordinates,
                                                   double x0, double y0,                    
                                                   size_t nX, size_t nY, 
                                                   double meanR, double sigmaR,
                                                   double meanShift, double sigmaShift,
                                                   size_t nL)
{
    size_t numAgents = generateHexPack(coordinates, x0, y0, nX, nY, meanR, nL);

    for (size_t iAgent=0; iAgent<numAgents; iAgent++) {
        
        const double shift = util::Random::getNextGauss(meanShift, sigmaShift);
        const double theta = util::Random::getNext(2.0*math::pi);
        
        coordinates[iAgent].pos[0] += shift*cos(theta);     // x
        coordinates[iAgent].pos[1] += shift*sin(theta);     // y
        coordinates[iAgent].resetPos0();
        coordinates[iAgent].r = util::Random::getNextGauss(meanR, sigmaR); // r
    }
    return numAgents;
}


void GeometryGenerator2D::mirror(AgentGeoVector &  coordinates, size_t idx)
{
    assert(idx<2); 
    
    const size_t numAgents = coordinates.size();
    for (size_t iAgent=0; iAgent<numAgents; iAgent++) {
        coordinates[iAgent].pos[idx] *= (-1.0);
        coordinates[iAgent].resetPos0();
    }
}

} // end namespace absmc
