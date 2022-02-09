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

#ifndef ABSMC_CELL_BASE_2D_H
#define ABSMC_CELL_BASE_2D_H

#include "core/agentBase.h"

namespace absmc {
   
class CellBase2D : public AgentBase<2> {
public:
    typedef Point<2, double> point_t;
public:    
    CellBase2D(AgentTypeId typeId, point_t pos, point_t  pos0, double r, point_t mobility)
        : AgentBase<2>(typeId, pos, pos0, r, mobility) { }
    
    virtual CellBase2D* clone() const = 0;  
    
    virtual double getScalarByName(std::string const& quantity) const {
        if (quantity=="wssOsi") return wssOsi;
        else if (quantity=="wssMax") return wssMax;
        else if (quantity=="drugConc") return drugConc;
        else return AgentBase<2>::getScalarByName(quantity);
    }
    
private:
    double          wssOsi, wssMax;
    double          drugConc;   
};
    

} // end namespace absmc

#endif
