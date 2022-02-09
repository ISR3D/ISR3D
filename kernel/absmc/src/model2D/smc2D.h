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

#ifndef ABSMC_SMC_2D_H
#define ABSMC_SMC_2D_H

#include "model2D/cellBase2D.h"

namespace absmc {

class SMC2D : public CellBase2D {
public:
    typedef AgentBase<2> base_t;
    typedef Point<2, double> point_t;
    enum CycleState { G0, G1, SG2M, A };
public:
     
    SMC2D(point_t pos, point_t pos0, double r, CycleState state_, size_t age_=0, size_t clock_=0) 
        : CellBase2D(tSMC2D, pos, pos0, r, point_t(1.0, 1.0) ), state(state_), age(age_), clock(clock_) { }
    
    virtual SMC2D* clone() const { return new SMC2D(*this); }
    
    void setAgentRule() { }
    virtual void execAgentRule(size_t iAgent, std::vector<base_t*> & agentsCreated, std::vector<base_t*> & agentsDeleted) { 
        // method not implemented yet
    }
    
    virtual double getValueByName(std::string const& quantity) const {
        if (quantity=="state") return static_cast<double>(state);
        else if (quantity=="age") return static_cast<double>(age);
        else if (quantity=="clock") return static_cast<double>(clock);
        else return CellBase2D::getScalarByName(quantity);
    }
    
private:
    CycleState      state;          // cell cycle state
    size_t          age;            // number of mitotic cycles completed
    size_t          clock;          // cell cycle clock
};
    

} // end namespace absmc

#endif
