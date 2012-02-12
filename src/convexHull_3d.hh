//----------------------------------------------------------------------------//
// 3D implementation of the convex hull algorithm.
// Based on the description by Timothy Chan:
//
// A minimalist's implementation of the 3-d divide-and-conquer convex hull algorithm
// http://www.cs.uwaterloo.ca/~tmchan
//----------------------------------------------------------------------------//
#ifndef __polytope_convexHull_3d__
#define __polytope_convexHull_3d__

namespace polytope {

namespace convexHull_helpers {
//------------------------------------------------------------------------------
// Timothy Chan    "ch3d.cc"    12/02    3-d lower hull (in C++)
//
// a simple implementation of the O(n log n) divide-and-conquer algorithm
//   input: coordinates of points
//     n x_0 y_0 z_0 ... x_{n-1} y_{n-1} z_{n-1}

//   output: indices of facets
//     i_1 j_1 k_1 i_2 j_2 k_2 ...

//   warning: ignores degeneracies and robustness
//   space: uses 6n pointers
//------------------------------------------------------------------------------
template<typename CoordType>
struct Point {
  CoordType x, y, z;
  Point *prev, *next;

  // These static members need to be initialized in a cc file.
  static CoordType INF;
  static Point     nil;
  static Point*    NIL;

  Point(): x(0), y(0), z(0), prev(0), next(0) {}
  Point(const CoordType xi, const CoordType yi, const CoordType zi,
        Point* previ, Point* nexti):
    x(xi), y(yi), z(zi), prev(previ), next(nexti) {}

  void act() {
    if (prev->next != this) {
      prev->next = next->prev = this;       // insert
    } else {
      prev->next = next; next->prev = prev; // delete
    }
  }
  
};

template<typename CoordType>
inline int turn(Point<CoordType> *p, Point<CoordType> *q, Point<CoordType> *r) {  // <0 iff cw
  if (p == Point<CoordType>::NIL || 
      q == Point<CoordType>::NIL ||
      r == Point<CoordType>::NIL) return 1;
  const double test = 
    (double(q->x) - double(p->x))*(double(r->y) - double(p->y)) -
    (double(r->x) - double(p->x))*(double(q->y) - double(p->y));
  return (test < 0.0 ? -1 : 1);
}

template<typename CoordType>
inline CoordType ch3dtime(Point<CoordType> *p, Point<CoordType> *q, Point<CoordType> *r) {  
  if (p == Point<CoordType>::NIL || 
      q == Point<CoordType>::NIL ||
      r == Point<CoordType>::NIL) return Point<CoordType>::INF;
  const CoordType den = turn(p, q, r);
  return (q->x - p->x)/den*(r->z - p->z) - (r->x - p->x)/den*(q->z - p->z);
}

template<typename CoordType>
Point<CoordType> *sort(Point<CoordType> P[], int n) {  // mergesort

  Point<CoordType> *a, *b, *c, head;

  if (n == 1) { P[0].next = Point<CoordType>::NIL; return P; }
  a = sort(P, n/2);
  b = sort(P+n/2, n-n/2); 
 c = &head;
  do
    if (a->x < b->x) { c = c->next = a; a = a->next; }
    else { c = c->next = b; b = b->next; }
  while (c != Point<CoordType>::NIL);
  return head.next;
}

template<typename CoordType>
void lowerHull(Point<CoordType> *list, 
               int n, 
               Point<CoordType> **A, 
               Point<CoordType> **B) {

  const CoordType INF = Point<CoordType>::INF;

  Point<CoordType> *u, *v, *mid;  
  CoordType t[6], oldt, newt;  
  int i, j, k, l, minl;

  if (n == 1) { A[0] = list->prev = list->next = Point<CoordType>::NIL; return; }

  for (u = list, i = 0; i < n/2-1; u = u->next, i++) ;
  mid = v = u->next;
  lowerHull(list, n/2, B, A);  // recurse on left and right sides
  lowerHull(mid, n-n/2, B+n/2*2, A+n/2*2);

  for ( ; ; )  // find initial bridge
    if (turn(u, v, v->next) < 0) v = v->next;
    else if (turn(u->prev, u, v) < 0) u = u->prev;  
    else break;

  // merge by tracking bridge uv over ch3dtime
  for (i = k = 0, j = n/2*2, oldt = -INF; ; oldt = newt) {  
    t[0] = ch3dtime(B[i]->prev, B[i], B[i]->next);  
    t[1] = ch3dtime(B[j]->prev, B[j], B[j]->next);    
    t[2] = ch3dtime(u, u->next, v);  
    t[3] = ch3dtime(u->prev, u, v);
    t[4] = ch3dtime(u, v->prev, v); 
    t[5] = ch3dtime(u, v, v->next);
    for (newt = INF, l = 0; l < 6; l++) 
      if (t[l] > oldt && t[l] < newt) { minl = l; newt = t[l]; }
    if (newt == INF) break;
    switch (minl) {
    case 0:  if (B[i]->x < u->x) A[k++] = B[i];  B[i++]->act();  break;
    case 1:  if (B[j]->x > v->x) A[k++] = B[j];  B[j++]->act();  break;
    case 2:  A[k++] = u = u->next;  break;
    case 3:  A[k++] = u;  u = u->prev;  break;
    case 4:  A[k++] = v = v->prev;  break;
    case 5:  A[k++] = v;  v = v->next;  break;
    }
  }
  A[k] = Point<CoordType>::NIL;

  u->next = v;  v->prev = u;  // now go back in ch3dtime to update pointers
  for (k--; k >= 0; k--) 
    if (A[k]->x <= u->x || A[k]->x >= v->x) {
      A[k]->act();
      if (A[k] == u) u = u->prev; else if (A[k] == v) v = v->next;
    }
    else { 
      u->next = A[k]; A[k]->prev = u; v->prev = A[k]; A[k]->next = v;
      if (A[k]->x < mid->x) u = A[k]; else v = A[k];
    }
}

}
}

#endif
