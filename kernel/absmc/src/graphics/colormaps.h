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

#ifndef ABSMC_COLORMAPS_H
#define ABSMC_COLORMAPS_H

#include <cassert>
#include <string>

#include "graphics/rgb.h"

namespace absmc {

namespace graphics {

class ColorMap {
public:
    virtual ~ColorMap() { }
    virtual rgb get(double x) const = 0;
};

class DefaultColorMap : public ColorMap {
public:
    DefaultColorMap() {
        cm[0]=rgb(0,0,0.529412);            cm[1]=rgb(0,0,0.560784);
        cm[2]=rgb(0,0,0.607843);            cm[3]=rgb(0,0,0.639216);
        cm[4]=rgb(0,0,0.686275);            cm[5]=rgb(0,0,0.717647);
        cm[6]=rgb(0,0,0.764706);            cm[7]=rgb(0,0,0.796078);
        cm[8]=rgb(0,0,0.843137);            cm[9]=rgb(0,0,0.87451);
        cm[10]=rgb(0,0,0.921569);           cm[11]=rgb(0,0,0.952941);
        cm[12]=rgb(0,0.00392157,0.994771);  cm[13]=rgb(0,0.027451,1);
        cm[14]=rgb(0,0.0745098,1);          cm[15]=rgb(0,0.105882,1);
        cm[16]=rgb(0,0.152941,1);           cm[17]=rgb(0,0.184314,1);
        cm[18]=rgb(0,0.231373,1);           cm[19]=rgb(0,0.262745,1);
        cm[20]=rgb(0,0.309804,1);           cm[21]=rgb(0,0.341176,1);
        cm[22]=rgb(0,0.388235,1);           cm[23]=rgb(0,0.419608,1);
        cm[24]=rgb(0,0.466667,1);           cm[25]=rgb(0,0.498039,1);
        cm[26]=rgb(0,0.545098,1);           cm[27]=rgb(0,0.576471,1);
        cm[28]=rgb(0,0.623529,1);           cm[29]=rgb(0,0.654902,1);
        cm[30]=rgb(0,0.701961,1);           cm[31]=rgb(0,0.733333,1);
        cm[32]=rgb(0,0.780392,1);           cm[33]=rgb(0,0.811765,1);
        cm[34]=rgb(0,0.858824,1);           cm[35]=rgb(0,0.890196,1);
        cm[36]=rgb(0,0.937255,1);           cm[37]=rgb(0,0.968627,1);
        cm[38]=rgb(0.0130719,1,0.994771);   cm[39]=rgb(0.0431373,1,0.968627);
        cm[40]=rgb(0.0901961,1,0.921569);   cm[41]=rgb(0.121569,1,0.890196);
        cm[42]=rgb(0.168627,1,0.843137);    cm[43]=rgb(0.2,1,0.811765);
        cm[44]=rgb(0.247059,1,0.764706);    cm[45]=rgb(0.278431,1,0.733333);
        cm[46]=rgb(0.32549,1,0.686275);     cm[47]=rgb(0.356863,1,0.654902);
        cm[48]=rgb(0.403922,1,0.607843);    cm[49]=rgb(0.435294,1,0.576471);
        cm[50]=rgb(0.482353,1,0.529412);    cm[51]=rgb(0.513725,1,0.498039);
        cm[52]=rgb(0.560784,1,0.45098);     cm[53]=rgb(0.592157,1,0.419608);
        cm[54]=rgb(0.639216,1,0.372549);    cm[55]=rgb(0.670588,1,0.341176);
        cm[56]=rgb(0.717647,1,0.294118);    cm[57]=rgb(0.74902,1,0.262745);
        cm[58]=rgb(0.796078,1,0.215686);    cm[59]=rgb(0.827451,1,0.184314);
        cm[60]=rgb(0.87451,1,0.137255);     cm[61]=rgb(0.905882,1,0.105882);
        cm[62]=rgb(0.952941,1,0.0588235);   cm[63]=rgb(0.984314,1,0.027451);
        cm[64]=rgb(1,0.984314,0);           cm[65]=rgb(1,0.952941,0);
        cm[66]=rgb(1,0.905882,0);           cm[67]=rgb(1,0.87451,0);
        cm[68]=rgb(1,0.827451,0);           cm[69]=rgb(1,0.796078,0);
        cm[70]=rgb(1,0.74902,0);            cm[71]=rgb(1,0.717647,0);
        cm[72]=rgb(1,0.670588,0);           cm[73]=rgb(1,0.639216,0);
        cm[74]=rgb(1,0.592157,0);           cm[75]=rgb(1,0.560784,0);
        cm[76]=rgb(1,0.513725,0);           cm[77]=rgb(1,0.482353,0);
        cm[78]=rgb(1,0.435294,0);           cm[79]=rgb(1,0.403922,0);
        cm[80]=rgb(1,0.356863,0);           cm[81]=rgb(1,0.32549,0);
        cm[82]=rgb(1,0.278431,0);           cm[83]=rgb(1,0.247059,0);
        cm[84]=rgb(1,0.2,0);                cm[85]=rgb(1,0.168627,0);
        cm[86]=rgb(1,0.121569,0);           cm[87]=rgb(1,0.0901961,0);
        cm[88]=rgb(1,0.0431373,0);          cm[89]=rgb(1,0.0130719,0);
        cm[90]=rgb(0.968627,0,0);           cm[91]=rgb(0.937255,0,0);
        cm[92]=rgb(0.890196,0,0);           cm[93]=rgb(0.858824,0,0);
        cm[94]=rgb(0.811765,0,0);           cm[95]=rgb(0.780392,0,0);
        cm[96]=rgb(0.733333,0,0);           cm[97]=rgb(0.701961,0,0);
        cm[98]=rgb(0.654902,0,0);           cm[99]=rgb(0.623529,0,0);
    }
    
    virtual rgb get(double x) const {
        assert(0<=x && x<=1);
        return cm[int(x*99)];
    }
private:
    rgb cm[100];
};

class BlackAndWhiteColorMap : public ColorMap {
public:
    virtual rgb get(double x) const {
        assert(0<=x && x<=1);
        return rgb(x, x, x);
    }
};

ColorMap* generateColorMap(std::string const& colorMapName) {
    if (colorMapName == "blackAndWhite") {
        return new BlackAndWhiteColorMap;
    } 
    else {
        return new DefaultColorMap;
    }
}


}  // namespace graphics

}  // namespace absmc

#endif
