#ifndef KDTREE_VISITOR_H
#define KDTREE_VISITOR_H

#include <cstddef>
#include <iostream>
#include "node.h"

namespace kdtree {

template<typename node_t>
class Visitor {
public:
    virtual ~Visitor() { }
    virtual void visit(node_t & node) = 0;
};


template<typename node_t>
class PrintVisitor : public Visitor<node_t> {
public:
    PrintVisitor(std::ostream & ostr_) : ostr(ostr_) { }

    virtual void visit(node_t & node) {
        ostr << "pIdxLo=" << node.pIdxLo << " pIdxHi=" << node.pIdxHi
             << " parent=" << node.parent << " bBox=" << node.bBox << std::endl;
    }
private:
    std::ostream & ostr;
};

} /* namespace kdtree */

#endif
