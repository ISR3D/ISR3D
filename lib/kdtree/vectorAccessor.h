#ifndef KDTREE_VECTOR_ACCESSOR_H
#define KDTREE_VECTOR_ACCESSOR_H

#include <vector>
#include "assert.h"
#include "geometry.h"

namespace kdtree {

/// \brief{Accessor template for use with KDTree, suitable to access objects stored by value in random access containers.}
/// VectorAccessor works with any random access container holding spatial objects by value,
/// i.e. it assumes that container_t defines operator[] to access the objects contained,
/// and object_t defines operator[] to access coordinates. This works fine e.g. with
/// container_type==std::vector<object_t>, object_type=kdtree::Point<k, float_t>.
/// To support different syntax for element and coordinate access, the user must provide his own accessor class.

template<size_t k, typename float_t, class container_t>
class VectorAccessor {
public:
    typedef Point<k, float_t> point_t;

    /// constructor
    VectorAccessor(container_t & data_) : data(data_) { }
    // compiler-generated copy constructor and assignment operator are fine.

    /// return number of objects stored in the underlying container
    size_t size() const {
        return data.size();
    }

    /// get coordinate iDim of element iData.
    float_t getCoordinate(size_t iData, size_t iDim) const {
        assert(iData < size());
        assert(iDim < k);
        return (data[iData])[iDim];
    }

    /// set coordinate iDim of element iData.
    void setCoordinate(size_t iData, size_t iDim, float_t value) {
        assert(iData < size());
        assert(iDim < k);
        (data[iData])[iDim] = value;
    }

    /// construct a Point data structure from spatial data accessible by acc.
    point_t makePoint(size_t iData) const {
        assert(iData < size());
        point_t p;
        for (size_t iDim=0; iDim<k; iDim++) {
            p.x[iDim] = getCoordinate(iData, iDim);
        }
        return p;
    }

private:
    container_t & data;
};


/// \brief{Accessor template for use with KDTree, suitable to access objects stored by pointer in random access containers.}

template<size_t k, typename float_t, class container_t>
class VectorPtrAccessor {
public:
    typedef Point<k, float_t> point_t;

    /// constructor
    VectorPtrAccessor(container_t & data_) : data(data_) { }
    // compiler-generated copy constructor and assignment operator are fine.

    /// return number of objects stored in the underlying container
    size_t size() const {
        return data.size();
    }

    /// get coordinate iDim of element iData.
    float_t getCoordinate(size_t iData, size_t iDim) const {
        assert(iData < size());
        assert(iDim < k);
        return (data[iData])->operator[](iDim);
    }

    /// set coordinate iDim of element iData.
    void setCoordinate(size_t iData, size_t iDim, float_t value) {
        assert(iData < size());
        assert(iDim < k);
        (data[iData])->operator[](iDim) = value;
    }

    /// construct a Point data structure from spatial data accessible by acc.
    point_t makePoint(size_t iData) const {
        assert(iData < size());
        point_t p;
        for (size_t iDim=0; iDim<k; iDim++) {
            p.x[iDim] = getCoordinate(iData, iDim);
        }
        return p;
    }

private:
    container_t & data;
};


} /* namespace kdtree */

#endif
