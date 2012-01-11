#ifndef POINT_HH
#define POINT_HH

#include <cmath>
#include "charybdisconfig.hh"
#include <algorithm>
#include <iostream>
#include <stdlib.h>

namespace Charybdis
{

//! \class Point
//! This class stores points in N-dimensional space.
template <int N>
class Point
{
  public:

  //! Initialize a zero point.
  Point() { std::fill(m_x, m_x + N, 0.0); }

  //! Initialize a point with all components set to the given value.
  explicit Point(Real val) { std::fill(m_x, m_x + N, val); }

  //! Initialize a point in two dimensions.
  Point(Real x1, Real x2) { m_x[0] = x1, m_x[1] = x2; }

  //! Initialize a point in three dimensions.
  Point(Real x1, Real x2, Real x3) { m_x[0] = x1, m_x[1] = x2, m_x[2] = x3; }

  //! Copy constructor.
  Point(const Point& x) { std::copy(x.m_x, x.m_x + N, m_x); }

  //! Destructor.
  ~Point() {}

  //! Assignment operator.
  Point& operator=(const Point& x)
  {
    if (&x != this)
    {
      std::copy(x.m_x, x.m_x + N, m_x);
    }
    return *this;
  }

  //! Creates and returns a random point within the bounding box
  //! whose limits are given by \a low and \a high.
  //! On UNIX, this function uses the random() function.
  static Point random(const Point& low, const Point& high)
  {
    Point x;
    for (int i = 0; i < N; ++i)
      x[i] = (double(::random())/RAND_MAX) * (high[i]-low[i]) + low[i];
    return x;
  }

  //! Creates and returns a random point within the hypercube [low,high]**N.
  //! On UNIX, this function uses the random() function.
  static Point random(Real low, Real high)
  {
    Point x;
    for (int i = 0; i < N; ++i)
      x[i] = (double(::random())/RAND_MAX) * (high-low) + low;
    return x;
  }

  //! Creates and returns a random point within the hypercube [0,1]**N.
  //! On UNIX, this function uses the random() function.
  static Point random()
  {
    Point x;
    for (int i = 0; i < N; ++i)
      x[i] = (double(::random())/RAND_MAX);
    return x;
  }

  //! Returns the corresponding component.
  Real operator[](int i) const { return m_x[i]; }
  Real& operator[](int i) { return m_x[i]; }

  //! Output operator.
  friend std::ostream& operator<<(std::ostream& s, const Point<N>& p)
  {
    s << "(";
    for (int i = 0; i < N; ++i)
    {
      s << p[i];
      if (i < (N-1))
        s << ", ";
    }
    s << ")";
    return s;
  }

  private:

  Real m_x[N];

};

//! Euclidean distance function.
template <int N>
Real distance(const Point<N>& p1, const Point<N>& p2)
{
  Real sum = 0.0;
  for (int i = 0; i < N; ++i)
    sum += (p2[i]-p1[i])*(p2[i]-p1[i]);
  return std::sqrt(sum);
}

}

#endif
