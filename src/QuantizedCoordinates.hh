#ifndef POLYTOPE_QUANTIZEDCOORDINATES_HH
#define POLYTOPE_QUANTIZEDCOORDINATES_HH

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "ReducedPLC.hh"
#include "Point.hh"
#include "polytope_geometric_utilities.hh"

namespace polytope {


//! \class QuantizedCoordinates - A Piecewise Linear Complex in 3D, or a 
//! Planar Straight Line Graph (PSLG) in 2D.
//! The quantized coordinates is a reduced PLC which
//!  a) stores lower and upper bounds for point quantization
//!  b) stores the center and radius of an inscribed bounding sphere
//!  c) has routines for passing points to/from floating point from/to integers
template<int Dimension, typename RealType>
class QuantizedCoordinates: public ReducedPLC<Dimension, RealType> {
public:

  typedef int64_t CoordHash;
  typedef Point2<CoordHash> IntPoint;
  typedef Point2<RealType> RealPoint;

  //! The bounding box size
  std::vector<RealType> low, high;

  //! The infinite bounding sphere parameters
  std::vector<RealType> center;
  RealType rinf;
   
  //! The quantized mesh spacing and degeneracy parameters
  RealType delta, mDegeneracy;
  CoordHash mCoordMax;

  //! Flag for checking if coordinates have been modified
  bool mCoordinatesModified;

  //------------------------------------------------------------------------
  //! Constructors
  //------------------------------------------------------------------------

  //! Default constructor
  QuantizedCoordinates():
    ReducedPLC<Dimension, RealType>(),
    low(),
    high(),
    center(){}
  QuantizedCoordinates(const QuantizedCoordinates& coords);

  //------------------------------------------------------------------------
  //! Initialize the coordinate system
  //------------------------------------------------------------------------
  void initialize(const std::vector<RealType>& allPoints) {
    // Initialize degeneracy and coordMax
    mDegeneracy = 1.5e-8;
    mCoordMax   = (1LL << 26);
    low.clear();  high.clear();  center.clear();
    
    // compute bounds for the box
    RealType plow[Dimension], phigh[Dimension];
    geometry::computeBoundingBox<Dimension,RealType>(allPoints, true, plow, phigh);
    
    // The infinite-radius bounding sphere
    RealType box[Dimension], cen[Dimension];
    rinf = 0.0;
    for (unsigned j = 0; j != Dimension; ++j) {
      box[j] =      phigh[j] - plow[j];
      cen[j] = 0.5*(phigh[j] + plow[j]);
      rinf   = std::max(rinf, 4.0*box[j]);
    }
    
    // Resize bounds to circumscribe sphere
    for (unsigned j = 0; j != Dimension; ++j) {
      plow [j] = cen[j] - rinf;
      phigh[j] = cen[j] + rinf;
    }

    // Finalize box parameters
    std::copy(plow , plow +Dimension, std::back_inserter(low   ));
    std::copy(phigh, phigh+Dimension, std::back_inserter(high  ));
    std::copy(cen  , cen  +Dimension, std::back_inserter(center));
    delta = std::max(mDegeneracy, 2.0*rinf/RealType(mCoordMax));

    // Fill in the ReducedPLC information
    std::vector<RealType> pts;
    std::vector<std::vector<int> > fcts;
    if (Dimension == 2) {
      pts.push_back(low [0]);  pts.push_back(low [1]);
      pts.push_back(high[0]);  pts.push_back(low [1]);
      pts.push_back(high[0]);  pts.push_back(high[1]);
      pts.push_back(low [0]);  pts.push_back(high[1]);

      fcts.resize(4, std::vector<int>(2));
      for (unsigned j = 0; j != 4; ++j) {
        fcts[j][0] = j;  fcts[j][1] = (j+1)%4;
      }
    } else {
      POLY_ASSERT(Dimension == 3);
      pts.push_back(low [0]);  pts.push_back(low [1]);  pts.push_back(low [2]);
      pts.push_back(high[0]);  pts.push_back(low [1]);  pts.push_back(low [2]);
      pts.push_back(high[0]);  pts.push_back(high[1]);  pts.push_back(low [2]);
      pts.push_back(low [0]);  pts.push_back(high[1]);  pts.push_back(low [2]);
      pts.push_back(low [0]);  pts.push_back(low [1]);  pts.push_back(high[2]);
      pts.push_back(high[0]);  pts.push_back(low [1]);  pts.push_back(high[2]);
      pts.push_back(high[0]);  pts.push_back(high[1]);  pts.push_back(high[2]);
      pts.push_back(low [0]);  pts.push_back(high[1]);  pts.push_back(high[2]);

      fcts.resize(6, std::vector<int>(4));
      fcts[0][0] = 0;  fcts[0][1] = 1;  fcts[0][2] = 2;  fcts[0][3] = 3;  // -z face
      fcts[0][0] = 4;  fcts[0][1] = 5;  fcts[0][2] = 6;  fcts[0][3] = 7;  // +z face
      fcts[0][0] = 0;  fcts[0][1] = 1;  fcts[0][2] = 5;  fcts[0][3] = 4;  // -y face
      fcts[0][0] = 2;  fcts[0][1] = 6;  fcts[0][2] = 7;  fcts[0][3] = 3;  // +y face
      fcts[0][0] = 0;  fcts[0][1] = 4;  fcts[0][2] = 7;  fcts[0][3] = 3;  // -x face
      fcts[0][0] = 1;  fcts[0][1] = 5;  fcts[0][2] = 6;  fcts[0][3] = 2;  // +x face
    }
    this->points = pts;
    this->facets = fcts;

    mCoordinatesModified = false;
  }

  //------------------------------------------------------------------------
  //! Clear and empty operations
  //------------------------------------------------------------------------
  void clear() {
    low.clear();  high.clear();  center.clear();  rinf=0;  delta=0;
  }
  bool empty() const {
     return (low.empty() and high.empty() and center.empty() and
             rinf == 0 and delta == 0);
  }

  //------------------------------------------------------------------------
  //! Degeneracy accessors
  //------------------------------------------------------------------------
  RealType getDegeneracy() {
    return mDegeneracy;
  }
  void setDegeneracy(const RealType degeneracy) {
    mDegeneracy = degeneracy;
    mCoordinatesModified = true;
  }

  //------------------------------------------------------------------------
  //! coordMax accessors
  //------------------------------------------------------------------------
  RealType getCoordMax() {
    return mCoordMax;
  }
  void setCoordMax(const CoordHash coordMax) {
    mCoordMax = coordMax;
    mCoordinatesModified = true;
  }

  //------------------------------------------------------------------------
  //! Expand bounding box based on new bounds
  //------------------------------------------------------------------------
  void expand(RealType plow [Dimension],
              RealType phigh[Dimension]) {
    // Pre-conditions
    POLY_ASSERT(!low.empty() and !high.empty() and !center.empty());
    POLY_ASSERT(!this->points.empty() and !this->facets.empty());

    // Resize the box
    RealType boxSize = 0.0;
    for (unsigned j = 0; j != Dimension; ++j) {
      plow [j] = std::min(low [j], plow [j]);
      phigh[j] = std::max(high[j], phigh[j]);
      boxSize  = std::max(boxSize, phigh[j]-plow[j]);
    }
    
    // Recompute a box center
    for (unsigned j = 0; j != Dimension; ++j) {
      center[j] = 0.5*(plow[j] + phigh[j]);
    }
    
    // Make a larger bounding sphere and resize the box to inscribe it
    rinf = 4.0*boxSize;
    for (unsigned j = 0; j != Dimension; ++j) {
      low [j] = std::min(low [j], center[j]-rinf);
      high[j] = std::max(high[j], center[j]+rinf);
    }
    
    // Recompute the quantized mesh spacing
    delta = std::max(mDegeneracy, 2.0*rinf/RealType(mCoordMax));
    mCoordinatesModified = false;
    
    // Update the bounding box vertices for the reduced PLC
    std::vector<RealType> pts;
    if (Dimension == 2) {
      pts.push_back(low [0]);  pts.push_back(low [1]);
      pts.push_back(high[0]);  pts.push_back(low [1]);
      pts.push_back(high[0]);  pts.push_back(high[1]);
      pts.push_back(low [0]);  pts.push_back(high[1]);
    } else {
      POLY_ASSERT(Dimension == 3);
      pts.push_back(low [0]);  pts.push_back(low [1]);  pts.push_back(low [2]);
      pts.push_back(high[0]);  pts.push_back(low [1]);  pts.push_back(low [2]);
      pts.push_back(high[0]);  pts.push_back(high[1]);  pts.push_back(low [2]);
      pts.push_back(low [0]);  pts.push_back(high[1]);  pts.push_back(low [2]);
      pts.push_back(low [0]);  pts.push_back(low [1]);  pts.push_back(high[2]);
      pts.push_back(high[0]);  pts.push_back(low [1]);  pts.push_back(high[2]);
      pts.push_back(high[0]);  pts.push_back(high[1]);  pts.push_back(high[2]);
      pts.push_back(low [0]);  pts.push_back(high[1]);  pts.push_back(high[2]);
    }
    this->points = pts;
  }

  //------------------------------------------------------------------------
  //! Expand bounding box based on new points
  //------------------------------------------------------------------------
  void expand(const std::vector<RealType>& newPoints,
              const bool globalReduce) {
    RealType plow[Dimension], phigh[Dimension];
    geometry::computeBoundingBox<Dimension,RealType>(newPoints, globalReduce, plow, phigh);
    expand(plow, phigh);
  }

  //------------------------------------------------------------------------
  //! Print out the bounding box info
  //------------------------------------------------------------------------
  friend std::ostream& operator<<(std::ostream& os, const QuantizedCoordinates& coords) {
    os << "Bounding Box  = ";
    for (unsigned j = 0; j != Dimension; ++j) {
      os << "(" << coords.low[j] << "," << coords.high[j] << ") ";
    }
    os << std::endl;
    os << "Sphere Center = (";
    for (unsigned j = 0; j != Dimension; ++j) {
       os << coords.center[j] << ",";
    }
    os << std::endl;
    os << "Sphere Radius = " << coords.rinf << std::endl;
    os << "Int Spacing   = " << coords.delta << std::endl;
    return os;
  }

  //------------------------------------------------------------------------
  //! Quantize a floating-point-precision point
  //------------------------------------------------------------------------
  inline
  IntPoint quantize(const RealType* pointIn) const {
    POLY_ASSERT(!this->mCoordinatesModified);
    return IntPoint(pointIn[0], pointIn[1], this->low[0], this->low[1], this->delta);
  }
   
  //------------------------------------------------------------------------
  //! Dequantize an integer point to floating-point-precision 
  //------------------------------------------------------------------------
  inline
  RealPoint dequantize(const IntPoint pointIn) const {
    POLY_ASSERT(!this->mCoordinatesModified);
    RealType x = (pointIn.x == mCoordMax) ? high[0] : pointIn.realx(low[0], delta);
    RealType y = (pointIn.y == mCoordMax) ? high[1] : pointIn.realy(low[1], delta);
    return RealPoint(x,y);
  }

  //------------------------------------------------------------------------
  //! Project a point to the bounding circle
  //------------------------------------------------------------------------
  inline
  RealPoint projectPoint(const RealType* rayPoint,
                         const RealType* rayDirection) const {
    RealPoint result;
    bool test = geometry::rayCircleIntersection(rayPoint, rayDirection, 
                                                &this->center[0], this->rinf, 
                                                this->mDegeneracy, &result.x);
    POLY_ASSERT(test);
    return result;
  }
};


// //! Partial specialization for 2D
// template <typename RealType>
// class QuantizedCoordinates<2, RealType>
// {
// public:

//   typedef int64_t CoordHash;
//   typedef Point2<CoordHash> IntPoint;
//   typedef Point2<RealType> RealPoint;

//   QuantizedCoordinates():
//     ReducedPLC<2, RealType>() {}

//   QuantizedCoordinates(const QuantizedCoordinates& coords);

//   void initialize(const std::vector<RealType> allPoints) {
//     this->initialize(allPoints);
//   }

//   //------------------------------------------------------------------------
//   //! Quantize a floating-point-precision point
//   //------------------------------------------------------------------------
//   inline
//   IntPoint quantize(const RealType* pointIn) {
//     POLY_ASSERT(!this->mCoordinatesModified);
//     return IntPoint(pointIn[0], pointIn[1], this->low[0], this->low[1], this->delta);
//   }
   
//   //------------------------------------------------------------------------
//   //! Dequantize an integer point to floating-point-precision 
//   //------------------------------------------------------------------------
//   inline
//   RealPoint dequantize(const IntPoint pointIn) {
//     POLY_ASSERT(!this->mCoordinatesModified);
//     return RealPoint(pointIn.realx(this->low[0], this->delta), 
//                      pointIn.realy(this->low[1], this->delta));
//   }

//   //------------------------------------------------------------------------
//   //! Project a point to the bounding circle
//   //------------------------------------------------------------------------
//   inline
//   RealPoint projectPoint(const RealPoint& rayPoint,
//                          const RealPoint& rayDirection) {
//     RealPoint result;
//     bool test = geometry::rayCircleIntersection(&rayPoint.x, &rayDirection.x, 
//                                                 &this->center[0], this->rinf, 
//                                                 this->mDegeneracy, &result.x);
//     POLY_ASSERT(test);
//     return result;
//   }
// };


// //! Partial specialization for 3D
// template <typename RealType>
// class QuantizedCoordinates<3, RealType>
// {
// public:

//   typedef int64_t CoordHash;
//   typedef Point3<CoordHash> IntPoint;
//   typedef Point3<RealType> RealPoint;

//   //------------------------------------------------------------------------
//   //! Quantize a floating-point-precision point
//   //------------------------------------------------------------------------
//   inline
//   IntPoint quantize(const RealType* pointIn) {
//     POLY_ASSERT(!this->mCoordinatesModified);
//     return IntPoint(pointIn[0], pointIn[1], pointIn[2],
//                     this->low[0], this->low[1], this->low[2], this->delta);
//   }
   
//   //------------------------------------------------------------------------
//   //! Dequantize an integer point to floating-point-precision 
//   //------------------------------------------------------------------------
//   inline
//   RealPoint dequantize(const IntPoint& pointIn) {
//     POLY_ASSERT(!this->mCoordinatesModified);
//     return RealPoint(pointIn.realx(this->low[0], this->delta),
//                      pointIn.realy(this->low[1], this->delta),
//                      pointIn.realz(this->low[2], this->delta));
//   }

//   //------------------------------------------------------------------------
//   //! Project a point to the bounding circle
//   //------------------------------------------------------------------------
//   inline
//   RealPoint projectPoint(const RealPoint& rayPoint,
//                          const RealPoint& rayDirection) {
//     RealPoint result;
//     bool test = geometry::raySphereIntersection(&rayPoint.x, &rayDirection.x, 
//                                                 &this->center[0], this->rinf, 
//                                                 this->mDegeneracy, &result.x);
//     POLY_ASSERT(test);
//     return result;
//   }
// };



} //end polytope namespace

#endif
