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
 
#ifndef ABSMC_PIXMAP_H
#define ABSMC_PIXMAP_H

#include <cstddef>
#include <cassert>

#include "graphics/rgb.h"

namespace absmc {

namespace graphics {

template<typename T>
class PixMap {
public:
    typedef T pixel_t;
public:
    PixMap(size_t nI_, size_t nJ_) : nI(nI_), nJ(nJ_), data(0) {
        assert(nI>0);
        assert(nJ>0);
        data = new pixel_t[nI*nJ];
    }
    
    ~PixMap() {
        delete[] data;
    }
    
    void setPixel(size_t i, size_t j, pixel_t const& pixel) {
        assert(0<=i && i<nI);
        assert(0<=j && i<nJ);
        data[i*nJ+j] = pixel;
    }
    
    pixel_t getPixel(size_t i, size_t j) const {
        assert(0<=i && i<nI);
        assert(0<=j && i<nJ);
        return data[i*nJ+j];
    }
    
    void setAll(pixel_t const& pixel) {
        for (size_t i=0; i<nI; i++) {
            for (size_t j=0; j<nJ; j++) {
                setPixel(i, j, pixel);
            }
        }
    }
    
private:
    PixMap(PixMap const& rhs);
    PixMap& operator=(PixMap const& rhs);
private:
    const size_t nI, nJ;
    pixel_t* data;
};

} // end namespace graphics

} // end namespace absmc

#endif
