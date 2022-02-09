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

#ifndef ABSMC_IEL_2D_H
#define ABSMC_IEL_2D_H

#include "model2D/cellBase2D.h"

namespace absmc {

class IEL2D : public CellBase2D {
public:
    typedef AgentBase<2> base_t;
    typedef Point<2, double> point_t;
public:
     
    IEL2D(point_t pos, point_t pos0, double r)
        : CellBase2D(tIEL2D, pos, pos0, r, point_t(1.0, 1.0) )  { }
    
    virtual IEL2D* clone() const { return new IEL2D(*this); }
    
    void setAgentRule() { }
    virtual void execAgentRule(size_t iAgent, std::vector<base_t*> & agentsCreated, std::vector<base_t*> & agentsDeleted) { 
        // method not implemented yet
    }
    
    virtual double getScalarByName(std::string const& quantity) const {
        return CellBase2D::getScalarByName(quantity);
    }
    
};
    

} // end namespace absmc

#endif
