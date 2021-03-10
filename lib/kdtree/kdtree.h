#ifndef KDTREE_KDTREE_H
#define KDTREE_KDTREE_H

#include <vector>
#include <limits>
#include <stack>
#include <queue>
#include <iostream>

#include "geometry.h"
#include "euclideanMetrics.h"
#include "node.h"
#include "visitor.h"
#include "hit.h"

namespace kdtree {

/// \brief{A kd-tree class for fast spatial searches on k-dimensional data.}
/// Read access to spatial coordinates of objects in the container is provided via the methods of the Accessor class.
/// The KDTree class does not modify the order of data it is constructed on top of. The data is not owned by the
/// KDTree; to ensure validity of the tree, data MUST NOT be modified as long as the KDTree object is in scope.

/// The KDTree uses a functor class named Accessor to access spatial coordinates of objects stored in the tree.
/// It can therefore be completely agnostic of how the user represents and stores these objects.
/// The Accessor class needs to be copyable; it must furthermore provide typedefs and define the methods
/// size(), getCoordinate() and setCoordinate() following the scheme below.

template<size_t k, typename float_t, class Accessor,
template<size_t _k, typename _float_t> class Metrics = EuclideanMetrics>
class KDTree {
public:
    typedef Point<k, float_t> point_t;
    typedef Box<k, float_t> box_t;
    typedef Node<k, float_t> node_t;
    typedef Hit<float_t> hit_t;

    typedef std::vector<size_t> idx_vec_t;

    typedef Metrics<k, float_t> metrics_t;
    typedef Visitor<node_t> visitor_t;

    /// Construct kd-tree from data accessable by Accessor acc.
    KDTree(Accessor acc_, metrics_t metrics_=metrics_t(), size_t terminalNodeSize_=2);

    /// Visit all nodes in depth-first order.
    void exploreDepthFirst(visitor_t & visitor);
    /// Visit all nodes in breadth-first order.
    void exploreBreadthFirst(visitor_t & visitor);
    /// Print raw node data (for debugging).
    void printRawNodes(std::ostream & ostr) const;

    // If the set of m nearest hits of a point is not unique (because there are several hits of the
    // same maximal distance included in the result set), the findMNearest() routines below can not guarantee
    // which of those hits are included in the result set.
    /// Find m nearest hits of an arbitrary point given by its coordinates in k-dim space;
    /// the result includes the query point if it equals one of those used to construct the tree.
    void findMNearest(point_t const& p, size_t m, idx_vec_t & results) const;
    /// Find m nearest hits of one of the points used to construct the tree;
    /// the result does not include the query point.
    void findMNearest(size_t dIdx, size_t m, idx_vec_t & results) const;

    /// Find all points within a ball around an arbitrary point given by its coordinates in k-dim space;
    /// the result includes the query point if it equals one of those used to construct the tree.
    void findWithinBall(point_t const& p, float_t radius, idx_vec_t & results) const;
    /// Find all points within a ball around one of the points used to construct the tree;
    /// the result does not include the query point.
    void findWithinBall(size_t dIdx, float_t radius, idx_vec_t & results) const;

    /// Find all points within a rectangular box in k-dim space.
    void findWithinBox(box_t const& box, idx_vec_t & results) const;

    // The convenience routines below are simply implemented in terms of the findWithinX() routines above,
    // so no performance benefit is achieved.
    /// Count number of points within a ball around an arbitrary point given by its coordinates in k-dim space.
    size_t countWithinBall(point_t const& p, float_t radius) const {
        idx_vec_t results;
        findWithinBall(p, radius, results);
        return static_cast<size_t>(results.size() );
    }
    /// Count number of points within a ball around one of the points used to construct the tree.
    size_t countWithinBall(size_t dIdx, float_t radius) const {
        idx_vec_t results;
        findWithinBall(dIdx, radius, results);
        return static_cast<size_t>(results.size() );
    }
    /// Count number of points within a rectangular box in k-dim space.
    size_t countWithinBox(box_t const& box) const {
        idx_vec_t results;
        findWithinBox(box, results);
        return static_cast<size_t>(results.size() );
    }
private:
    KDTree(KDTree const& rhs);
    KDTree& operator=(KDTree const& rhs);

    static const size_t I_NO_CHILD;
    static const size_t I_NO_IGNORE;

    // Return index of node containing an arbitrary point given by its coordinates in k-dim space.
    // Note: This does currently not work correctly for periodic coordinate space.
    size_t locate(point_t const& p) const;
    // Return index of node containing one of the points used to construct the tree.
    size_t locate(size_t dIdx) const;

    float_t initSearchList(point_t const& p, size_t pIdxToIgnore, size_t m,
                            std::priority_queue<hit_t> & pqueue, size_t initNode) const;
    void implFindMNearest(point_t const& p, size_t pIdxToIgnore, size_t m,
                          idx_vec_t & results, size_t initNode) const;
    void implFindWithinBall(point_t const& p, size_t pIdxToIgnore, float_t radius, idx_vec_t & results) const;

private:
    Accessor acc;                     // functor to access underlying data
    metrics_t metrics;                // class defining metrics.

    const size_t nPts;                // number of points
    const size_t terminalNodeSize;    // capacity of terminal nodes
    size_t nNodesEst, nNodes;         // estimated and constructed number of nodes

    // terminology: xIndex refers to a mapping(type vector), xIdx refers to a key (type size_t).
    // dIndex = direct index into data; data[dIdx] = dIdx'th point as seen by the user.
    // pIndex = index used internally for partitioning; every node stores a range pIdxLo<=pIdx<=pIdxHi.
    //          data[dIdx=pIndex[pIdx]] = data corresponding to a point stored at the node.
    // rIndex = reverse index;
    //          data[pIdx=rIndex[dIdx]] = dIdx'th point as seen by the user.
    std::vector<size_t> pIndex;
    std::vector<size_t> rIndex;

    std::vector<node_t> nodes;         // the tree
};

} /* namespace kdtree */

#endif
