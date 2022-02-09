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
 
#ifndef ABSMC_GEOMETRY_GENERATOR_2D_H
#define ABSMC_GEOMETRY_GENERATOR_2D_H

#include <vector>

#include "core/geometry.h"

namespace absmc {

struct AgentGeometry2D {
    
    typedef Point<2, double> point_t;
    
    point_t pos, pos0;
    double r;
    
    AgentGeometry2D(point_t pos_, point_t pos0_, double r_) : pos(pos_), pos0(pos0), r(r_) { }
    
    AgentGeometry2D(double x, double y, double x0, double y0, double r_) : r(r_) {
        pos[0]=x; pos[1]=y;
        pos0[0]=x0; pos0[1]=y0;
    }
    
    void resetPos0() { pos0 = pos; }
}; 
    
/// Generate a regular geometry with hexagonal structure.
class GeometryGenerator2D {
public:
    typedef std::vector<AgentGeometry2D> AgentGeoVector;
public:
    /// nX is the number of equilibrium-spaced cells that fit into lX
    /// such that the outer cells have equilibrium overlap with the boundaries.
    static double getLx(size_t nI, double R) {
        const double s  = sqrt(3.0)*R;
        const double lX = nI * s;
        return lX;
    }
    /// nY is the number of equilibrium-spaced cells that fit into lY
    /// such that the outer cells touch the boundaries in one point.
    static double getLy(size_t nJ, double R) {
        const double lY = (1.5*nJ + 0.5) * R ;
        return lY;    
    }
    /// Generate a regular structure of nX*nY cells with radius R and origin in (x0,y0). 
    /// nL=2..5 is the number of small lining layer cells per large cell.
    /// Output is an IC-list style vector of coordinates (x1,y1,z1,r1,x2,y2,...).
    static size_t generateHexPack(AgentGeoVector & coordinates,
                                double x0, double y0,                    
                                size_t nX, size_t nY, double R,
                                size_t nL);
    
    /// Generate hexpack-based structure, but with Gaussian-distributed radius and shift from hex-grid locations.
    static size_t generatePerturbedHexPack(AgentGeoVector &  coordinates,
                                         double x0, double y0,                    
                                         size_t nX, size_t nY, 
                                         double meanR, double sigmaR,
                                         double meanShift, double sigmaShift,
                                         size_t nL);
    
    /// Given a regular hexpack-based structure, generate a lining layer of nL small cells per large cell.
    /// The lower-most layer of large cells is assumed to be stored in the first nX data sets of cellCoordinates.
    static size_t generateLiningLayer(AgentGeoVector &  liningCoordinates,
                                    AgentGeoVector &  cellCoordinates,
                                    size_t nX, size_t nL);
    
    
    /// Flip coordinate indexed by idx.
    static void mirror(AgentGeoVector & coordinates, size_t idx);
    
};

} // end namespace absmc

#endif
