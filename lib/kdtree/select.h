#ifndef KDTREE_SELECT_H
#define KDTREE_SELECT_H

#include <cstddef>
#include <vector>

namespace kdtree {

/// \brief Selection of m-th smallest value by partitioning, done via index.
/// Permutes index[iLo..iHi] such that
/// data_iDim[index[iLo..iLo+m-1]] <= data_iDim[index[iLo+m]] <= data_iDim[index[iLo+m+1..iHi]]
/// and returns index[m]. Data in the container is not permuted.
/// data_iDim[j] denotes the iDim-th spatial coordinate of the j-th data element in the container.
template<typename float_t, class Accessor>
size_t selectViaIndex(size_t m, size_t iLo, size_t iHi,  std::vector<size_t> & index, Accessor acc, size_t iDim)
{
    assert(iLo <= iHi);

    // number of elements in the range to be partitioned
    const size_t N = iHi-iLo+1;
    // shifted index; working on sIndex instead of index saves adding offset iLo at every operation
    size_t* sIndex = &index[iLo];   //

    size_t curLo = 0;
    size_t curHi = N-1;

    while (true) {
        // if current partition contains <= 2 elements: done
        if (curHi <= curLo+1) {
            if ( (curHi==curLo+1) &&
                    (acc.getCoordinate(sIndex[curHi], iDim) < acc.getCoordinate(sIndex[curLo], iDim)) ) {
                std::swap(sIndex[curLo], sIndex[curHi]);
            }
            return sIndex[m];
        }
        // otherwise work on the active partition.
        else {
            // choose median of left, right and center elements as partitioning element, indexed by sIndex[curLo+1],
            // and make sure that the "sentinel condition" acc(curLo) <= acc(curLo+1) <= acc(curHi) holds.
            const size_t curMed = (curLo+curHi) >> 1;
            std::swap(sIndex[curMed], sIndex[curLo+1]);
            if (acc.getCoordinate(sIndex[curLo], iDim) > acc.getCoordinate(sIndex[curHi], iDim) )
                { std::swap(sIndex[curLo], sIndex[curHi]); }
            if (acc.getCoordinate(sIndex[curLo+1], iDim) > acc.getCoordinate(sIndex[curHi], iDim) )
                { std::swap(sIndex[curLo+1], sIndex[curHi]); }
            if (acc.getCoordinate(sIndex[curLo], iDim) > acc.getCoordinate(sIndex[curLo+1], iDim) )
                { std::swap(sIndex[curLo], sIndex[curLo+1]); }

            // partitioning element
            const size_t iA = sIndex[curLo+1];
            const float_t a = acc.getCoordinate(iA, iDim);
            // initialise pointers for scanning up and down
            size_t i = curLo+1;
            size_t j = curHi;

            while (true) {
                do i++; while (acc.getCoordinate(sIndex[i], iDim) < a); // scan up
                do j--; while (acc.getCoordinate(sIndex[j], iDim) > a); // scan down
                if (j < i) break;
                std::swap(sIndex[i], sIndex[j]);
            }
            // move partitioning element to final position
            sIndex[curLo+1] = sIndex[j];
            sIndex[j] = iA;
            // keep partitioning containing the m-th element active
            if (j >= m) curHi = j-1;
            if (j <= m) curLo = i;
        }

    } // end of main loop
}


} /* namespace kdtree */

#endif
