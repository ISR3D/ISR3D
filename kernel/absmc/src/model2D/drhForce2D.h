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

#ifndef ABSMC_DRH_FORCE_2D_H
#define ABSMC_DRH_FORCE_2D_H

#include "core/geometry.h"
#include "core/euclideanMetrics.h"
#include "model2D/cellBase2D.h"
#include "util/mathUtil.h"

namespace absmc {
    
struct DRHForceConstants {
    static const double E;
    static const double nu;
    static const double E1;
    static const double K;
    static const double piTimesE;
    static const double piTimesE1Over4;
};
   
// numerical constants used by DRH cell-cell force, cell-boundary force and hoop force:
const double DRHForceConstants::E  = 1.0;                              // Young's modulus
const double DRHForceConstants::nu = 0.0;                              // Poisson's ratio
const double DRHForceConstants::E1 = E / (1.0-nu*nu);                  // effective Young's modulus for plane strain
const double DRHForceConstants::K  = math::pi/8.0*E1;                  // attraction coefficent;
const double DRHForceConstants::piTimesE = math::pi*E;                 // helper
const double DRHForceConstants::piTimesE1Over4 = math::pi*E1/4.0;      // helper   
   
/// Hoop force.   
class UnaryForceDRH {
public:
    typedef AgentBase<2> agent_t;
    typedef Point<2, double> point_t;
    
    UnaryForceDRH(double crFactor_) : crFactor(crFactor_) { }
    
    void accumulate(agent_t & agent) { 
    
        const double y = agent.getPos()[1];
        double y0 = agent.getPos0()[1];
        
        // work-around for very small y0: inner torus radius cannot be smaller than its thickness
        if (fabs(y0)<agent.getR() ) { y0 = math::sgn(y0)*agent.getR(); }
    
        // linear elastic force
        // new by 09/03/09: configurable relative hoop stiffness in compressive regime (1 = tensile regime)
        double factor = (fabs(y)-fabs(y0) > 0.0) ? 1.0 : crFactor;
        const double k = factor * DRHForceConstants::piTimesE * agent.getRSqr() / (y0*y0);
        const double fR = k*(y0-y);
        
        #ifdef ABSMC_EXCESSIVE_CHECKING
        if (isnan(fR)) { std::cerr << "isnan(fR)\n"; }
        #endif
            
        agent.addForce(point_t(0.0, fR));
    }
    
private:
    double crFactor;
};

/// Cell-cell force based on Hooke's model of elastic cylinders.
class BinaryForceDRH {
public:
    typedef AgentBase<2> agent_t;
    typedef Point<2, double> point_t;
    typedef EuclideanMetrics<2, double> metrics_t;
    
    BinaryForceDRH() { 
        // it is fine to initialise EuclideanMetrics with its default constructor
    }
    
    void accumulate(agent_t & agent1, agent_t const& agent2) {
    
        // check if cells are in contact
        
        const double sSqr = metrics.distSqr(agent1.getPos(), agent2.getPos());
        const double s = sqrt(sSqr);
        const double R1 = agent1.getR();
        const double R2 = agent2.getR();
        if (s > R1+R2 ) {
//             force[0] = force[1] = 0.0;
            return;
        }
        
        // compute "contact radius"
        const double R1Sqr = R1*R1;
        const double R2Sqr = R2*R2;
        const double sMaxSqr = fabs(R1Sqr-R2Sqr); // distance at which contact radius is maximal
        double aSqr = 0.0; // contact radius
        
        if (sSqr < sMaxSqr+1e-6) { // regularization introduced by BS 16/09/2008
            // note: epsilon term preliminarily added 14/01/09 
            // to work around instability at R1-R2->0, s->0
             aSqr = std::min(R1Sqr, R2Sqr);
        }
        else { // original term by DRH
            aSqr = ( R1Sqr*R2Sqr - 0.25 * (sSqr-R1Sqr-R2Sqr)*(sSqr-R1Sqr-R2Sqr) ) / sSqr;
        }

        // Hertzian repulsive force
        double f12 = - DRHForceConstants::piTimesE1Over4 * aSqr * (1.0/R1 + 1.0/R2);
        // attractive force
        f12 += 2.0 * DRHForceConstants::K * sqrt(aSqr);
        
        #ifdef ABSMC_EXCESSIVE_CHECKING
        if (isnan(f12)) { std::cerr << "drh_force: isnan(f12) : \n"; }
        if (isnan(1.0/s)) std::cerr << "morse_force: isnan(1/s)\n";
        #endif
        
        point_t const& pos1 = agent1.getPos();
        point_t const& pos2 = agent2.getPos();
        
//         force[0] = (pos2[0]-pos1[0]) / s * f12;
//         force[1] = (pos2[1]-pos1[1]) / s * f12;
        
        agent1.addForce(point_t((pos2[0]-pos1[0]) / s * f12, (pos2[1]-pos1[1]) / s * f12) );
    }
private:
    metrics_t metrics; 
};
    

} // end namespace absmc

#endif
