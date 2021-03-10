#ifndef KDTREE_KDTREE_HH
#define KDTREE_KDTREE_HH

#include "kdtree.h"
#include "select.h"

#include <algorithm>
#include <functional>
#include <limits>
#include <cmath>
#include <cassert>

using std::cout;
using std::endl;

namespace {

/// Task struct used to implement recursion via a heap of pending tasks.
struct Task {
    size_t nodeIdx;
    size_t cutDim;
    Task(size_t nodeIdx_, size_t cutDim_) : nodeIdx(nodeIdx_), cutDim(cutDim_) { }
};

/// Function object to compare indices by some coordinate of the data they refer to.
template<class Accessor>
class CmpIndexByReferencedData : public std::binary_function<size_t, size_t, bool> {
public:
    CmpIndexByReferencedData(Accessor acc_, size_t iDim_) : acc(acc_), iDim(iDim_) { }
    bool operator()(size_t dIdx1, size_t dIdx2) {
        return (acc.getCoordinate(dIdx1, iDim) < acc.getCoordinate(dIdx2, iDim) );
    }
private:
    Accessor acc;
    size_t iDim;
};

} /* end unnamed namespace */


namespace kdtree{

// definition of static data members
template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
const size_t
KDTree<k, float_t, Accessor, Metrics>::I_NO_CHILD = 0;

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
const size_t
KDTree<k, float_t, Accessor, Metrics>::I_NO_IGNORE = std::numeric_limits<size_t>::max();


template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
KDTree<k, float_t, Accessor, Metrics>::KDTree(Accessor acc_, metrics_t metrics_, size_t terminalNodeSize_)
    : acc(acc_),
      metrics(metrics_),
      nPts(acc.size()),
      terminalNodeSize(terminalNodeSize_),
      nNodesEst(0), nNodes(0),
      pIndex(nPts), rIndex(nPts)
{
    assert(nPts > 0);
    assert(terminalNodeSize > 0);

    assert(nPts < I_NO_IGNORE);

    // estimate number of nodes
    if (nPts <= terminalNodeSize) { // trivial case
        nNodesEst = 1;
    }
    else if (terminalNodeSize==2){ // we know the exact number of nodes a priori
        size_t m = 1;
        for (size_t tmp=nPts ; tmp; tmp >>=1) { m <<= 1; }
        nNodesEst = std::min(m-1, 2*nPts-m/2-1);
    }
    else { // otherwise take a safe estimate
        // number of terminal nodes required to store nPts points
        const size_t nTReq = static_cast<size_t>(ceil(static_cast<float_t>(nPts) / terminalNodeSize) );
        // smallest power of 2 not less than nTReq
        size_t nTFul = 1;
        for (size_t tmp=nTReq-1 ; tmp; tmp >>=1) { nTFul <<= 1; }
        // total number of nodes in completely balanced tree where terminal level has nTFul nodes
        nNodesEst = 2*nTFul-1;
    }

    nodes.reserve(nNodesEst);

    // initialise index
    for (size_t dIdx=0; dIdx<nPts; dIdx++) pIndex[dIdx]=dIdx;

    // create root node
    const float_t maxVal = std::numeric_limits<float_t>::max() - 1;
    point_t lo, hi;
    for (size_t iDim=0; iDim<k; iDim++) {
        lo[iDim] = -maxVal;
        hi[iDim] = +maxVal;
    }
    node_t root(0, nPts-1, 0, I_NO_CHILD, I_NO_CHILD, box_t(lo, hi) );
    // insert into tree and push task onto stack
    size_t iNode = 0;
    nodes[iNode] = root;

    typedef Task task_t;
    std::stack<task_t> tasks;

    if (nPts>terminalNodeSize) {
        tasks.push(task_t(0, 0) );
    }

    // main loop over pending tasks
    while (!tasks.empty() ) {
        task_t const& task = tasks.top();
        const size_t iCurNode = task.nodeIdx;
        const node_t& curNode = nodes[iCurNode];
        const size_t cutDim = task.cutDim;
        tasks.pop();

        // get range of pIdx of current node
        size_t pIdxLo = curNode.pIdxLo;
        size_t pIdxHi = curNode.pIdxHi;
        // compute number of points in node and offset of median point
        size_t nP = pIdxHi-pIdxLo+1;
        size_t pidxMed = (nP-1) >> 1;

        // actual partitioning:
        // this permutes pIndex[pIdxLo..pIdxHi] such that
        // data[pIndex[pIdxLo..pIdxLo+pidxMed-1]] <= data[pIndex[pIdxLo+pidxMed]] <= data[pIndex[pIdxLo+pidxMed+1..pIdxHi]]
        // and returns pIndex[k].
        selectViaIndex<float_t, Accessor>(pidxMed, pIdxLo, pIdxHi, pIndex, acc, cutDim);

        // the STL algorithm could possibly replace the home-cooked function...
//         typedef size_t* iterator_t;
//         typedef CmpIndexByReferencedData<accessor_t> cmp_t;
//         std::nth_element<iterator_t, cmp_t> (&pIndex[pIdxLo], &pIndex[pIdxMed], &pIndex[pIdxHi+1], cmp_t(acc, cutDim) );

        // compute bounding boxes for child nodes
        const float_t cutVal = acc.getCoordinate(pIndex[pIdxLo+pidxMed], cutDim);
        box_t leftBBox(curNode.bBox.lo, curNode.bBox.hi);
        leftBBox.hi.x[cutDim] = cutVal;
        box_t rightBBox(curNode.bBox.lo, curNode.bBox.hi);
        rightBBox.lo.x[cutDim] = cutVal;

        // create child nodes and insert into tree;
        // median point (data[pIndex[pIdxLo+pidxMed]]) goes into the left node.
        node_t leftNode = node_t(pIdxLo, pIdxLo+pidxMed, iCurNode, I_NO_CHILD, I_NO_CHILD, leftBBox);
        node_t rightNode = node_t(pIdxLo+pidxMed+1, pIdxHi, iCurNode, I_NO_CHILD, I_NO_CHILD, rightBBox);

        nodes[++iNode] = leftNode;
        nodes[iCurNode].left = iNode;
        nodes[++iNode] = rightNode;
        nodes[iCurNode].right = iNode;
        // if newly created nodes contain more than terminalNodeSize points, create and push tasks to split them.
        if (pidxMed+1>terminalNodeSize) {
            tasks.push(task_t(iNode-1, (cutDim+1)%k) );
        }
        if (nP-1-pidxMed>terminalNodeSize) {
            tasks.push(task_t(iNode, (cutDim+1)%k) );
        }
    } // end of main loop
    nNodes = iNode+1;

    // create reverse index
    for (size_t dIdx=0; dIdx<nPts; dIdx++) rIndex[pIndex[dIdx]]=dIdx;

}

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
inline float_t
KDTree<k, float_t, Accessor, Metrics>::initSearchList(
        point_t const& p, size_t pIdxToIgnore, size_t m, std::priority_queue<hit_t> & pqueue, size_t iNode) const
{
    // node used for initialisation must contain at least m points
    assert(nodes[iNode].pIdxHi-nodes[iNode].pIdxLo+1 >= m);

    // add first m points to search list unconditionally...
    size_t pIdx;
    for (pIdx=nodes[iNode].pIdxLo; static_cast<size_t>(pqueue.size())<m; pIdx++) {
        if (pIdx==pIdxToIgnore) continue;

        const float_t dSqr = metrics.distSqr(p, acc.makePoint(pIndex[pIdx]) );
        pqueue.push(hit_t(pIndex[pIdx], dSqr) );
    }
    float_t maxDistSqrSoFar = pqueue.top().distSqr;

    // now loop over remaining points, replacing the most distant point in the result list
    // with the currently examined one, if it is closer
    for ( ; pIdx<=nodes[iNode].pIdxHi; pIdx++) {
        if (pIdx==pIdxToIgnore) continue;

        const float_t dSqr = metrics.distSqr(p, acc.makePoint(pIndex[pIdx]) );
        // if q is closer to p than the most distant point in the result list,
        // replace the latter by q and update search radius.
        if (dSqr < maxDistSqrSoFar) {
            pqueue.pop();
            pqueue.push(hit_t(pIndex[pIdx], dSqr) );
            maxDistSqrSoFar = pqueue.top().distSqr;
        }
    }
    return maxDistSqrSoFar;
}

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
inline void
KDTree<k, float_t, Accessor, Metrics>::implFindMNearest(
        point_t const& p, size_t pIdxToIgnore, size_t m, idx_vec_t & results, size_t initNode) const
{
    // search list: a priority queue
    // sorted by distance from p such that the most distant point is on top.
    std::priority_queue<hit_t> searchList;

    // initialise search list with m closest points contained in that box.
    float_t maxDistSqrSoFar = initSearchList(p, pIdxToIgnore, m, searchList, initNode);
    float_t maxDistSoFar = sqrt(maxDistSqrSoFar);

    typedef Task task_t;
    std::stack<task_t> tasks;

    // now search the tree starting from root, like for the r-ball search,
    // but updating the search-radius whenever a point on the search list is replaced.
    // if the search list was initialized using points from the root node, nothing to do here.
    if (initNode != 0) {
        tasks.push(task_t(0, 0) );
    }

    while (!tasks.empty() ) {
        task_t const& task = tasks.top();
        const size_t iCurNode = task.nodeIdx;
        const node_t& curNode = nodes[iCurNode];
        const size_t cutDim = task.cutDim;
        tasks.pop();

        // if current node is a terminal node
        if (curNode.left == I_NO_CHILD) {
             // examine all points contained in terminal node
            for (size_t pIdx=curNode.pIdxLo; pIdx<=curNode.pIdxHi; pIdx++) {
                if (pIdx==pIdxToIgnore) continue;

                const float_t dSqr = metrics.distSqr(p, acc.makePoint(pIndex[pIdx]) );
                // if q is closer to p than the most distant point in the result list,
                // replace the latter by q and update search radius.
                if (dSqr < maxDistSqrSoFar) {
                    searchList.pop();
                    searchList.push(hit_t(pIndex[pIdx], dSqr) );
                    maxDistSqrSoFar = searchList.top().distSqr;
                    maxDistSoFar = sqrt(maxDistSqrSoFar);
                }
            }
        }
        // otherwise, add tasks for child nodes if they are closer than r.
        else {
            const size_t iLeftNode = curNode.left;
            const size_t iRightNode = curNode.right;

            // skip node used for initialisation; check only nodes that may contain closer points.
            box_t const& lBBox = nodes[iLeftNode].bBox;
            if ((iLeftNode!=initNode) &&
                 metrics.intersects(cutDim, lBBox.lo.x[cutDim], lBBox.hi.x[cutDim], p.x[cutDim]-maxDistSoFar, p.x[cutDim]+maxDistSoFar)) {
                tasks.push(task_t(iLeftNode, (cutDim+1)%k) );
            }
            box_t const& rBBox = nodes[iRightNode].bBox;
            if ((iRightNode!=initNode) &&
                 metrics.intersects(cutDim, rBBox.lo.x[cutDim], rBBox.hi.x[cutDim], p.x[cutDim]-maxDistSoFar, p.x[cutDim]+maxDistSoFar)) {
                tasks.push(task_t(iRightNode, (cutDim+1)%k) );
            }
        }
    }

    // inefficient: copy result from priority queue to vector; TODO: avoid.
    while(!searchList.empty() ) {
        results.push_back(searchList.top().dIdx );
        searchList.pop();
    }
}

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
void
KDTree<k, float_t, Accessor, Metrics>::findMNearest(
        point_t const& p, size_t m, idx_vec_t & results) const
{
    // useless to search for more neighbours than points in the tree
    assert(m <= nPts);

    results.clear();
    if (m==0) return;

    // locate box containing p, then traverse tree upwards
    // to find smallest box containing enough points to init search list.
    size_t initNode = locate(p);
    while (nodes[initNode ].pIdxHi-nodes[initNode ].pIdxLo+1 < m) {
        initNode = nodes[initNode].parent;
    }

    implFindMNearest(p, I_NO_IGNORE, m, results, initNode);
}

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
void
KDTree<k, float_t, Accessor, Metrics>::findMNearest(
        size_t dIdx, size_t m, idx_vec_t & results) const
{
    // useless to search for more neighbours than points in the tree, not counting the query point
    assert(m <= nPts-1);

    results.clear();
    if (m==0) return;

    const point_t p = acc.makePoint(dIdx);
    const size_t pIdxToIgnore = rIndex[dIdx];

    // locate box containing data[idx],
    // then traverse tree upwards to find smallest box containing enough points to initialise search list.
    size_t initNode = locate(dIdx);
    while (nodes[initNode].pIdxHi-nodes[initNode].pIdxLo+1 < m+1) {
        initNode = nodes[initNode].parent;
    }

    implFindMNearest(p, pIdxToIgnore, m, results, initNode);
}


template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
inline void
KDTree<k, float_t, Accessor, Metrics>::implFindWithinBall(
        point_t const& p, size_t pIdxToIgnore, float_t radius, idx_vec_t & results) const
{
    results.clear();

    // start searching at root node only if it is closer to p than r.
    const float_t rSqr = radius*radius;
    if (metrics.distSqr(nodes[0].bBox, p) > rSqr) return;

    typedef Task task_t;
    std::stack<task_t> tasks;
    tasks.push(task_t(0, 0) );

    while (!tasks.empty() ) {
        task_t const& task = tasks.top();
        const size_t iCurNode = task.nodeIdx;
        const node_t& curNode = nodes[iCurNode];
        const size_t cutDim = task.cutDim;
        tasks.pop();

        // if current node is a terminal node
        if (curNode.left == I_NO_CHILD) {
             // examine all points contained in terminal node
            for (size_t pIdx=curNode.pIdxLo; pIdx<=curNode.pIdxHi; pIdx++) {
                if (pIdx==pIdxToIgnore) continue;

                const float_t dSqr = metrics.distSqr(p, acc.makePoint(pIndex[pIdx]) );
                if (dSqr <= rSqr) {
                    results.push_back(pIndex[pIdx]);
                }
            }
        }
        // otherwise, add tasks for child nodes if they are closer than r.
        else {
            const size_t iLeftNode = curNode.left;
            const size_t iRightNode = curNode.right;

            // need to compare only the coordinate along which the current node is split.
            box_t const& lBBox = nodes[iLeftNode].bBox;
            if (metrics.intersects(cutDim, lBBox.lo.x[cutDim], lBBox.hi.x[cutDim], p.x[cutDim]-radius, p.x[cutDim]+radius)) {
                tasks.push(task_t(iLeftNode, (cutDim+1)%k) );
            }
            box_t const& rBBox = nodes[iRightNode].bBox;
            if (metrics.intersects(cutDim, rBBox.lo.x[cutDim], rBBox.hi.x[cutDim], p.x[cutDim]-radius, p.x[cutDim]+radius)) {
                tasks.push(task_t(iRightNode, (cutDim+1)%k) );
            }
        }
    }
}

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
void
KDTree<k, float_t, Accessor, Metrics>::findWithinBall(
        point_t const& p, float_t radius, idx_vec_t & results) const
{
    implFindWithinBall(p, I_NO_IGNORE, radius, results);
}

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
void
KDTree<k, float_t, Accessor, Metrics>::findWithinBall(
        size_t dIdx, float_t radius, idx_vec_t & results) const
{
    const point_t p = acc.makePoint(dIdx);
    const size_t pIdxToIgnore = rIndex[dIdx];

    implFindWithinBall(p, pIdxToIgnore, radius, results);
}

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
void
KDTree<k, float_t, Accessor, Metrics>::findWithinBox(
        box_t const& box, idx_vec_t & results) const
{
    results.clear();

    typedef Task task_t;
    std::stack<task_t> tasks;
    // query box is always contained in root_box, so we can start at root node.
    tasks.push(task_t(0, 0) );

    while (!tasks.empty() ) {
        task_t const& task = tasks.top();
        const size_t iCurNode = task.nodeIdx;
        const node_t& curNode = nodes[iCurNode];
        const size_t cutDim = task.cutDim;
        tasks.pop();

        // if current node is a terminal node
        if (curNode.left == I_NO_CHILD) {
            // examine all points contained in terminal node
            for (size_t pIdx=curNode.pIdxLo; pIdx<=curNode.pIdxHi; pIdx++) {
                const float_t dSqr = metrics.distSqr(box, acc.makePoint(pIndex[pIdx]) );
                if (dSqr==0) {
                    results.push_back(pIndex[pIdx]);
                }
            }
        }
        // otherwise, add tasks for child nodes if they are closer than r.
        else {
            const size_t iLeftNode = curNode.left;
            const size_t iRightNode = curNode.right;
            // need to compare only the coordinate along which the current node is split.
            // box overlaps with bounding box of left child node
            box_t const& lBBox = nodes[iLeftNode].bBox;
            if (metrics.intersects(cutDim, lBBox.lo[cutDim], lBBox.hi[cutDim], box.lo[cutDim], box.hi[cutDim]) ) {
                tasks.push(task_t(iLeftNode, (cutDim+1)%k) );
            }
            // box overlaps with bounding box of right child node
            box_t const& rBBox = nodes[iRightNode].bBox;
            if (metrics.intersects(cutDim, rBBox.lo[cutDim], rBBox.hi[cutDim], box.lo[cutDim], box.hi[cutDim]) ) {
                tasks.push(task_t(iRightNode, (cutDim+1)%k) );
            }
        }
    }
}


template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
void
KDTree<k, float_t, Accessor, Metrics>::exploreDepthFirst(visitor_t & visitor)
{
    size_t start = 0; // start exploring at root node

    std::stack<size_t> stak;
    stak.push(start);

    while(!stak.empty() ) {
        node_t & cur = nodes[stak.top()];
        stak.pop();
        visitor.visit(cur);
        if (cur.right != I_NO_CHILD) {
            stak.push(cur.right);
        }
        if (cur.left != I_NO_CHILD) {
            stak.push(cur.left);
        }
    }
}

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
void
KDTree<k, float_t, Accessor, Metrics>::exploreBreadthFirst(visitor_t & visitor)
{
    size_t start = 0; // start exploring at root node

    std::queue<size_t> next;
    next.push(start);

    while(!next.empty() ) {
        node_t & cur = nodes[next.front()];
        next.pop();
        visitor.visit(cur);
        if (cur.left != I_NO_CHILD) {
            next.push(cur.left);
        }
        if (cur.right != I_NO_CHILD) {
            next.push(cur.right);
        }
    }
}

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
void
KDTree<k, float_t, Accessor, Metrics>::printRawNodes(std::ostream & ostr) const
{
    for (size_t iNode=0; iNode<nNodes; iNode++) {
        node_t const& node = nodes[iNode];
        ostr << "iNode=" << iNode << "  parent=" << node.parent
             << "  left=" << node.left  << "  right=" << node.right
             << "  pIdxLo=" << node.pIdxLo    << "  pIdxHi=" << node.pIdxHi
             << "  bBox=" << node.bBox;
        if (node.left == I_NO_CHILD) {
            ostr << "  terminal node: ";
            for (size_t pIdx=node.pIdxLo; pIdx<=node.pIdxHi; pIdx++) {
                ostr << pIndex[pIdx] << " ";
            }
        }
        ostr << std::endl;
    }
}

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
inline size_t
KDTree<k, float_t, Accessor, Metrics>::locate(point_t const& p) const
{
    size_t iNode = 0; // start at root
    size_t iDim = 0;

    while (true) {
        const size_t iLeft = nodes[iNode].left;
        if (iLeft == I_NO_CHILD) return iNode;

        if (p.x[iDim] <= nodes[iLeft].bBox.hi.x[iDim]) {
            iNode = iLeft;
        }
        else {
            iNode = nodes[iNode].right;
        }
        iDim = (iDim+1) % k;
    }
}

template<size_t k, typename float_t, class Accessor, template<size_t _k, typename _float_t> class Metrics>
inline size_t
KDTree<k, float_t, Accessor, Metrics>::locate(size_t dIdx) const
{
    const size_t pIdx = rIndex[dIdx];

    size_t iNode = 0; // start at root

    while (true) {
        const size_t iLeft = nodes[iNode].left;
        if (iLeft == I_NO_CHILD) return iNode;

        if (pIdx <= nodes[iLeft].pIdxHi) {
            iNode = iLeft;
        } else {
            iNode = nodes[iNode].right;
        }
    }
}

} /* namespace kdtree */

#endif
