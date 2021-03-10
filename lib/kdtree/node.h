#ifndef KDTREE_NODE_H
#define KDTREE_NODE_H

#include "geometry.h"

namespace kdtree {

/// Data structure representing a node of the kd-tree.
template<size_t k, typename float_t>
struct Node {
    size_t pIdxLo, pIdxHi;
    size_t parent, left, right;
    Box<k, float_t> bBox;

    Node(size_t pIdxLo_, size_t pIdxHi_, size_t parent_, size_t left_, size_t right_, Box<k, float_t> const& bBox_)
        : pIdxLo(pIdxLo_), pIdxHi(pIdxHi_), parent(parent_), left(left_), right(right_), bBox(bBox_) { }
};

} /* namespace kdtree */

#endif
