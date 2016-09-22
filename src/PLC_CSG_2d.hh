//------------------------------------------------------------------------------
// 2D adaptation of the 3D CSG methods in PLC_CSG_3d.hh.
//
// Thinly adapted from the following CSG libraries:
//
// Original CSG.JS library by Evan Wallace (http://madebyevan.com), under the MIT license.
// GitHub: https://github.com/evanw/csg.js/
// 
// C++ port by Tomasz Dabrowski (http://28byteslater.com), under the MIT license.
// GitHub: https://github.com/dabroz/csgjs-cpp/
//------------------------------------------------------------------------------
#ifndef __Polytope_PLC_CSG_2d__
#define __Polytope_PLC_CSG_2d__

#include <map>
#include <vector>
#include <algorithm>

#include "ReducedPLC.hh"
#include "Point.hh"
#include "polytope_geometric_utilities.hh"

namespace polytope {
namespace CSG {

//------------------------------------------------------------------------------
// Public interface methods (forward declaration).
//------------------------------------------------------------------------------
template<typename RealType> ReducedPLC<2, RealType> csg_union    (const ReducedPLC<2, RealType>& a, const ReducedPLC<2, RealType>& b);
template<typename RealType> ReducedPLC<2, RealType> csg_intersect(const ReducedPLC<2, RealType>& a, const ReducedPLC<2, RealType>& b);
template<typename RealType> ReducedPLC<2, RealType> csg_subtract (const ReducedPLC<2, RealType>& a, const ReducedPLC<2, RealType>& b);

//------------------------------------------------------------------------------
// Everything from here down is implementation detail.
//
// PointType vector methods.
//------------------------------------------------------------------------------
namespace CSG_internal_2d {
template<typename RealType> Point2<RealType> lerp(const Point2<RealType>& a, const Point2<RealType>& b, const double v) { 
  return Point2<RealType>(a.x + (b.x - a.x)*v,
                          a.y + (b.y - a.y)*v);
}
template<typename RealType> RealType         length2(const Point2<RealType>& a) { return a.x*a.x + a.y*a.y; }
template<typename RealType> RealType         dot(const Point2<RealType>& a, const Point2<RealType>& b) { return geometry::dot<2, RealType>(&a.x, &b.x); }
template<typename RealType> Point2<RealType> perp(const Point2<RealType>& a, const Point2<RealType>& b) { Point2<RealType> result; geometry::cross<2, RealType>(&a.x, &b.x, &result.x); return result; }

//------------------------------------------------------------------------------
// Represents a vertex of a segment. 
//------------------------------------------------------------------------------
template<typename RealType>
struct Vertex {
  typedef Point2<RealType> PointType;
  PointType pos;
  Vertex(): pos() {}
  Vertex(const PointType& posi): pos(posi) {}
};

//------------------------------------------------------------------------------
// Forward declarations.
//------------------------------------------------------------------------------
template<typename RealType> struct Line;
template<typename RealType> struct Segment;
template<typename RealType> struct Node;

//------------------------------------------------------------------------------
// Represents a line in 2D space.
//------------------------------------------------------------------------------
template<typename RealType>
struct Line {
  typedef Point2<RealType> PointType;
  PointType direction, normal;
  double w;
  static double EPSILON;

  Line(): direction(), normal(), w(0.0) {}
  Line(const PointType& a, const PointType& b): 
    direction(b - a),
    normal(direction.y, -(direction.x)),
    w(dot(normal, a)) {}
  bool ok() const { return (length2(this->normal) > 0.0) && (length2(this->direction) > 0.0); }
  void flip() { direction = -direction; normal = -normal; w *= -1.0; }
  void splitSegment(const Segment<RealType> & segment,
                    std::vector<Segment<RealType> > & coplanarFront,
                    std::vector<Segment<RealType> > & coplanarBack,
                    std::vector<Segment<RealType> > & front, 
                    std::vector<Segment<RealType> > & back) const;
};

//------------------------------------------------------------------------------
// Represents a line segment.
//------------------------------------------------------------------------------
template<typename RealType>
struct Segment {
  std::vector<Vertex<RealType> > vertices;
  Line<RealType> line;

  void flip() {
    std::swap(vertices[0], vertices[1]);
    line.flip();
  }

  Segment() {};
  Segment(const std::vector<Vertex<RealType> >& list):
    vertices(list), 
    line(list[0].pos, list[1].pos) {
    POLY_ASSERT(list.size() == 2);
  }
};

//------------------------------------------------------------------------------
// Holds a node in a BSP tree. A BSP tree is built from a collection of segments
// by picking a segment to split along. That segment (and all other colinear
// segments) are added directly to that node and the other segments are added to
// the front and/or back subtrees. This is not a leafy BSP tree since there is
// no distinction between internal and leaf nodes.
//------------------------------------------------------------------------------
template<typename RealType> 
struct Node {
  std::vector<Segment<RealType> > segments;
  Node<RealType> * front;
  Node<RealType> * back;
  Line<RealType> line;

  Node(): front(0), back(0) {}
  Node(const std::vector<Segment<RealType> > & list);
  ~Node() { delete front; delete back; }

  Node * clone() const;
  void clipTo(const Node * other);
  void invert();
  void build(const std::vector<Segment<RealType> > & segment);
  std::vector<Segment<RealType> > clipSegments(const std::vector<Segment<RealType> > & list) const;
  std::vector<Segment<RealType> > allSegments() const;
};

//------------------------------------------------------------------------------
// Vertex implementation
//------------------------------------------------------------------------------
// Invert all orientation-specific data. Called when the orientation of a segment is flipped.
template<typename RealType> Vertex<RealType> flip(const Vertex<RealType>& v) { return v; }

// Create a new vertex between this vertex and `other` by linearly
// interpolating all properties using a parameter of `t`. Subclasses should
// override this to interpolate additional properties.
template<typename RealType>
Vertex<RealType> interpolate(const Vertex<RealType>& a, const Vertex<RealType>& b, double t) {
  Vertex<RealType> ret;
  ret.pos = lerp(a.pos, b.pos, t);
  return ret;
}

//------------------------------------------------------------------------------
// Line implementation
//------------------------------------------------------------------------------
// Split `segment` by this line if needed, then put the segment or segment
// fragments in the appropriate lists.  Colinear segments go into either
// `collinear Front` or `collinearBack` depending on their orientation with
// respect to this line. Segments in front or in back of this line go into
// either `front` or `back`.
template<typename RealType>
void Line<RealType>::splitSegment(const Segment<RealType>& segment, 
                                   std::vector<Segment<RealType> >& collinearFront, 
                                   std::vector<Segment<RealType> >& collinearBack,
                                   std::vector<Segment<RealType> > & front, 
                                   std::vector<Segment<RealType> > & back) const {
  enum {
    COLLINEAR = 0,
    FRONT = 1,
    BACK = 2,
    SPANNING = 3
  };

  // Classify each point as well as the entire segment into one of the above
  // four classes.
  int segmentType = 0;
  std::vector<int> types;

  for (size_t i = 0; i < 2; i++) {
    double t = dot(this->normal, segment.vertices[i].pos) - this->w;
    int type = (t < -EPSILON) ? BACK : ((t > EPSILON) ? FRONT : COLLINEAR);
    segmentType |= type;
    types.push_back(type);
  }

  // Put the segment in the correct list, splitting it when necessary.
  switch (segmentType) {
  case COLLINEAR: 
    {
      if (dot(this->normal, segment.line.normal) > 0)
        collinearFront.push_back(segment);
      else 
        collinearBack.push_back(segment);
      break;
    }
  case FRONT: 
    {
      front.push_back(segment);
      break;
    }
  case BACK: 
    {
      back.push_back(segment);
      break;
    }
  case SPANNING: 
    {
      std::vector<Vertex<RealType> > f, b;
      int ti = types[0], tj = types[1];
      Vertex<RealType> vi = segment.vertices[0], vj = segment.vertices[1];
      if (ti != BACK) f.push_back(vi);
      if (ti != FRONT) b.push_back(vi);
      if (tj != BACK) f.push_back(vj);
      if (tj != FRONT) b.push_back(vj);
      if ((ti | tj) == SPANNING) {
        double t = (this->w - dot(this->normal, vi.pos)) / dot(this->normal, vj.pos - vi.pos);
        Vertex<RealType> v = interpolate(vi, vj, t);
        f.push_back(v);
        b.push_back(v);
      }
      if (f.size() == 2) front.push_back(Segment<RealType>(f));
      if (b.size() == 2) back.push_back(Segment<RealType>(b));
      break;
    }
  }
}

//------------------------------------------------------------------------------
// Node implementation
//------------------------------------------------------------------------------
// Return a new CSG solid representing space in either this solid or in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
template<typename RealType>
Node<RealType>* csg_union(const Node<RealType>* a1, const Node<RealType>* b1) {
  Node<RealType> * a = a1->clone();
  Node<RealType> * b = b1->clone();
  a->clipTo(b);
  b->clipTo(a);
  b->invert();
  b->clipTo(a);
  b->invert();
  a->build(b->allSegments());
  Node<RealType>* ret = new Node<RealType>(a->allSegments());
  delete a; a = 0;
  delete b; b = 0;
  return ret;
}

// Return a new CSG solid representing space in this solid but not in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
template<typename RealType>
Node<RealType>* csg_subtract(const Node<RealType>* a1, const Node<RealType> * b1) {
  Node<RealType>* a = a1->clone();
  Node<RealType>* b = b1->clone();
  a->invert();
  a->clipTo(b);
  b->clipTo(a);
  b->invert();
  b->clipTo(a);
  b->invert();
  a->build(b->allSegments());
  a->invert();
  Node<RealType>* ret = new Node<RealType>(a->allSegments());
  delete a; a = 0;
  delete b; b = 0;
  return ret;
}

// Return a new CSG solid representing space both this solid and in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
template<typename RealType>
Node<RealType>* csg_intersect(const Node<RealType>* a1, const Node<RealType>* b1) {
  Node<RealType>* a = a1->clone();
  Node<RealType>* b = b1->clone();
  a->invert();
  b->clipTo(a);
  b->invert();
  a->clipTo(b);
  b->clipTo(a);
  a->build(b->allSegments());
  a->invert();
  Node<RealType>* ret = new Node<RealType>(a->allSegments());
  delete a; a = 0;
  delete b; b = 0;
  return ret;
}

// Convert solid space to empty space and empty space to solid space.
template<typename RealType>
void Node<RealType>::invert() {	
  for (size_t i = 0; i < this->segments.size(); i++) this->segments[i].flip();
  this->line.flip();
  if (this->front) this->front->invert();
  if (this->back) this->back->invert();
  std::swap(this->front, this->back);
}

// Recursively remove all segments in `segments` that are inside this BSP
// tree.
template<typename RealType>
std::vector<Segment<RealType> > 
Node<RealType>::clipSegments(const std::vector<Segment<RealType> >& list) const {
  if (!this->line.ok()) return list;
  std::vector<Segment<RealType> > list_front, list_back;	
  for (size_t i = 0; i < list.size(); i++) {
    this->line.splitSegment(list[i], list_front, list_back, list_front, list_back);
  }
  if (this->front) list_front = this->front->clipSegments(list_front);
  if (this->back) list_back = this->back->clipSegments(list_back);
  else list_back.clear();
	
  list_front.insert(list_front.end(), list_back.begin(), list_back.end());
  return list_front;
}

// Remove all segments in this BSP tree that are inside the other BSP tree
// `bsp`.
template<typename RealType>
void Node<RealType>::clipTo(const Node<RealType>* other) {
  this->segments = other->clipSegments(this->segments);
  if (this->front) this->front->clipTo(other);
  if (this->back) this->back->clipTo(other);
}

// Return a list of all segments in this BSP tree.
template<typename RealType>
std::vector<Segment<RealType> > Node<RealType>::allSegments() const {
  std::vector<Segment<RealType> > list = this->segments;
  std::vector<Segment<RealType> > list_front, list_back;
  if (this->front) list_front = this->front->allSegments();
  if (this->back) list_back = this->back->allSegments();
  list.insert(list.end(), list_front.begin(), list_front.end());
  list.insert(list.end(), list_back.begin(), list_back.end());
  return list;
}

template<typename RealType>
Node<RealType>* Node<RealType>::clone() const {
  Node<RealType>* ret = new Node<RealType>();
  ret->segments = this->segments;
  ret->line = this->line;
  if (this->front) ret->front = this->front->clone();
  if (this->back) ret->back = this->back->clone();
  return ret;
}

// Build a BSP tree out of `segments`. When called on an existing tree, the
// new segments are filtered down to the bottom of the tree and become new
// nodes there. Each set of segments is partitioned using the first segment
// (no heuristic is used to pick a good split).
template<typename RealType>
void Node<RealType>::build(const std::vector<Segment<RealType> > & list) {
  if (!list.size()) return;
  if (!this->line.ok()) this->line = list[0].line;
  std::vector<Segment<RealType> > list_front, list_back;
  for (size_t i = 0; i < list.size(); i++) {
    this->line.splitSegment(list[i], this->segments, this->segments, list_front, list_back);
  }
  if (list_front.size()) {
    if (!this->front) this->front = new Node<RealType>;
    this->front->build(list_front);
  }
  if (list_back.size()) {
    if (!this->back) this->back = new Node<RealType>;
    this->back->build(list_back);
  }
}

template<typename RealType>
Node<RealType>::Node(const std::vector<Segment<RealType> > & list): front(0), back(0) {
  build(list);
}

//------------------------------------------------------------------------------
// Convert between polytope::ReducedPLC <-> vector<Segment>
//------------------------------------------------------------------------------
// Convert a ReducedPLC to a set of CSG segments.
template<typename RealType>
std::vector<Segment<RealType> >
ReducedPLCtoSegments(const ReducedPLC<2, RealType>& model) {
  typedef Point2<RealType> PointType;
  const unsigned nfacets = model.facets.size();
  const std::vector<RealType>& coords = model.points;
  std::vector<Segment<RealType> > list;
  std::vector<Vertex<RealType> > seg(2);
  for (size_t i = 0; i != nfacets; ++i) {
    POLY_ASSERT(model.facets[i].size() == 2);
    const unsigned k1 = model.facets[i][0], k2 = model.facets[i][1];
    seg[0] = Vertex<RealType>(PointType(coords[2*k1], coords[2*k1+1]));
    seg[1] = Vertex<RealType>(PointType(coords[2*k2], coords[2*k2+1]));
    list.push_back(Segment<RealType>(seg));
  }
  POLY_ASSERT(list.size() == model.facets.size());
  return list;
}

// Convert a set of segments to a ReducedPLC.
// Note we do not remove degeneracies here, so it's up to the caller
// to do with that as they will!
template<typename RealType>
ReducedPLC<2, RealType>
ReducedPLCfromSegments(const std::vector<Segment<RealType> >& segments) {
  // typedef Point2<RealType> PointType;
  ReducedPLC<2, RealType> result;
  const unsigned nfacets = segments.size();
  POLY_ASSERT(nfacets >= 3);
  for (size_t i = 0; i != nfacets; ++i) {
    POLY_ASSERT(segments[i].vertices.size() == 2);
    result.facets.push_back(std::vector<int>());
    result.facets.back().push_back(result.points.size()/2);
    result.points.push_back(segments[i].vertices[0].pos.x);
    result.points.push_back(segments[i].vertices[0].pos.y);
    result.facets.back().push_back(result.points.size()/2);
    result.points.push_back(segments[i].vertices[1].pos.x);
    result.points.push_back(segments[i].vertices[1].pos.y);
  }
  return result;
}

} // end namespace CSG_internal

//------------------------------------------------------------------------------
// Public interface implementation
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<2, RealType> csg_union(const ReducedPLC<2, RealType>& a,
                                  const ReducedPLC<2, RealType>& b) {
  using namespace CSG_internal_2d;
  Node<RealType>* A = new Node<RealType>(ReducedPLCtoSegments(a));
  Node<RealType>* B = new Node<RealType>(ReducedPLCtoSegments(b));
  Node<RealType>* AB = csg_union(A, B);                                  // <-- only difference
  std::vector<Segment<RealType> > segments = AB->allSegments();
  // delete A; A = 0;
  // delete B; B = 0;
  // delete AB; AB = 0;
  return ReducedPLCfromSegments(segments);
}

template<typename RealType>
ReducedPLC<2, RealType> csg_intersect(const ReducedPLC<2, RealType>& a,
                                      const ReducedPLC<2, RealType>& b) {
  using namespace CSG_internal_2d;
  Node<RealType>* A = new Node<RealType>(ReducedPLCtoSegments(a));
  Node<RealType>* B = new Node<RealType>(ReducedPLCtoSegments(b));
  Node<RealType>* AB = csg_intersect(A, B);                           // <-- only difference
  std::vector<Segment<RealType> > segments = AB->allSegments();
  delete A; A = 0;
  delete B; B = 0;
  delete AB; AB = 0;
  return ReducedPLCfromSegments(segments);
}

template<typename RealType>
ReducedPLC<2, RealType> csg_subtract(const ReducedPLC<2, RealType>& a,
                                     const ReducedPLC<2, RealType>& b) {
  using namespace CSG_internal_2d;
  Node<RealType>* A = new Node<RealType>(ReducedPLCtoSegments(a));
  Node<RealType>* B = new Node<RealType>(ReducedPLCtoSegments(b));
  Node<RealType>* AB = csg_subtract(A, B);                            // <-- only difference
  std::vector<Segment<RealType> > segments = AB->allSegments();
  delete A; A = 0;
  delete B; B = 0;
  delete AB; AB = 0;
  return ReducedPLCfromSegments(segments);
}

} // end namespace CGS
} // end namespace polytope

#endif
