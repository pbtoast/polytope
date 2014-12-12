#ifndef POLYTOPE_QUANTIZEDCOORDINATES_HH
#define POLYTOPE_QUANTIZEDCOORDINATES_HH

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "ReducedPLC.hh"
#include "Point.hh"
#include "DimensionTraits.hh"
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
  typedef typename DimensionTraits<Dimension, RealType>::IntPoint  IntPoint;
  typedef typename DimensionTraits<Dimension, RealType>::RealPoint RealPoint;
  
  //! The bounding box size
  RealPoint low, high;

  //! Infinite bounding sphere
  RealPoint center;
  RealType rinf;

  //! Parameters for controlling degeneracy spacing
  RealType mDegeneracy;
  CoordHash mCoordMax;

  //! The degeneracy spacing
  RealType delta;

  //! Flag for checking if coordinates have been modified
  bool mCoordinatesModified;
  bool mOuterBox;

  //------------------------------------------------------------------------
  //! Constructors
  //------------------------------------------------------------------------

  //! Default constructor
  QuantizedCoordinates():
    ReducedPLC<Dimension, RealType>(),
    low(), 
    high(), 
    center(),
    mCoordinatesModified(false),
    mOuterBox(false) {}
   
  //! Constuctor + initialize
  QuantizedCoordinates(const std::vector<RealType>& points):
    ReducedPLC<Dimension, RealType>(),
    low(), 
    high(),
    center(),
    mCoordinatesModified(false),
    mOuterBox(false)
  {
    initialize(points);
    POLY_ASSERT(rinf   != 0.0);
    POLY_ASSERT(delta  != 0.0);
    POLY_ASSERT(!mCoordinatesModified);
  }


  //------------------------------------------------------------------------
  //! Initialize the coordinate system
  //------------------------------------------------------------------------
  void initialize(const std::vector<RealType>& allPoints,
                  const RealType degeneracy) {
    
    // Initialize degeneracy and coordMax
    if (!mCoordinatesModified) {
      mDegeneracy = degeneracy;
      delta = degeneracy;
    }
    
    // compute bounds for the box
    geometry::computeBoundingBox<Dimension,RealType>(allPoints, true, &low.x, &high.x);

    // The infinite-radius bounding sphere
    RealPoint box;
    rinf = 0.0;
    for (unsigned j = 0; j != Dimension; ++j) {
      center[j] = 0.5*(high[j] + low[j]);
      box[j]    =      high[j] - low[j];
      rinf      = std::max(rinf, 0.75*box[j]);
      POLY_ASSERT(center[j] - rinf < low [j] and
                  center[j] + rinf > high[j]);
    }

    // Resize bounds to circumscribe sphere
    for (unsigned j = 0; j != Dimension; ++j) {
      low [j] = std::min(low [j], center[j] - rinf);
      high[j] = std::max(high[j], center[j] + rinf);
    }

    // Figure out the maximum integer coordinate needed
    mCoordMax = CoordHash(ceil(rinf / mDegeneracy));
    POLY_ASSERT2(mCoordMax < std::numeric_limits<CoordHash>::max(),
                 "Quantization Error: The floating point degeneracy tolerance " <<
                 "specified of " << mDegeneracy << " overflows the integer grid. " <<
                 "Please specify a coarser tolerance or use a longer integer.");
    
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
    rinf=0; delta=0;
    mOuterBox=false; mCoordinatesModified=false;
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

  // //------------------------------------------------------------------------
  // //! coordMax accessors
  // //------------------------------------------------------------------------
  // RealType getCoordMax() {
  //   return mCoordMax;
  // }
  // void setCoordMax(const CoordHash coordMax) {
  //   mCoordMax = coordMax;
  //   mCoordinatesModified = true;
  // }

  //------------------------------------------------------------------------
  //! Expand bounding box based on new bounds
  //------------------------------------------------------------------------
  void expand(const RealType* plow, const RealType* phigh) {
    // Pre-conditions
    POLY_ASSERT(!this->points.empty() and !this->facets.empty());

    // Resize the box
    RealType boxHalf = 0.0;
    for (unsigned j = 0; j != Dimension; ++j) {
      POLY_ASSERT(plow [j] <= phigh[j]);
      POLY_ASSERT(center[j] <= std::max(phigh[j], high[j]));
      POLY_ASSERT(center[j] >= std::min(plow [j], low [j]));
      boxHalf = std::max(boxHalf, std::max(phigh[j], high[j]) - center[j]);
      boxHalf = std::max(boxHalf, center[j] - std::min(plow [j], low [j]));
      low [j] = std::min(low [j], plow [j]);
      high[j] = std::max(high[j], phigh[j]);
    }
    
    // Circumscribe the new box dimensions with the infinite sphere
    rinf = std::max(rinf, 1.5*boxHalf);
    
    // Inflate box dimensions to inscribe the new infinite sphere 
    for (unsigned j = 0; j != Dimension; ++j) {
      low [j] = std::min(low [j], center[j]-rinf);
      high[j] = std::max(high[j], center[j]+rinf);
      POLY_ASSERT(low[j] <= high[j]);
    }

#if HAVE_MPI
    for (unsigned j = 0; j != Dimension; ++j) {
      low [j] = allReduce(low [j], MPI_MIN, MPI_COMM_WORLD);
      high[j] = allReduce(high[j], MPI_MAX, MPI_COMM_WORLD);
    }
    rinf = 0.5*(high[0] - low[0]);
#endif

    // Recompute the maximum coordinate
    mCoordMax = CoordHash(ceil(rinf / mDegeneracy));
    POLY_ASSERT2(mCoordMax < std::numeric_limits<CoordHash>::max(),
                 "Quantization Error: After expanding the bounding box to represent " <<
                 "Voronoi nodes, the floating point degeneracy tolerance " <<
                 mDegeneracy << " now overflows the integer grid. " <<
                 "Please specify a coarser tolerance or use a longer integer.");
    
    mCoordinatesModified = false;

    // Post conditions
    POLY_ASSERT(rinf > 0);
  }

  //------------------------------------------------------------------------
  //! Expand bounding box based on new points
  //------------------------------------------------------------------------
  void expand(const std::vector<RealType>& newPoints,
              const bool globalReduce) {
    RealType plow[Dimension], phigh[Dimension];
    geometry::computeBoundingBox<Dimension,RealType>(newPoints, 
						     globalReduce, 
						     plow, phigh);
    expand(plow, phigh);
  }

  //------------------------------------------------------------------------
  //! Print out the bounding box info
  //------------------------------------------------------------------------
  friend std::ostream& operator<<(std::ostream& os, const QuantizedCoordinates& coords) {
    os << "Bounding Box       = ";
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
  //! Check if a point is inside inner bounding box
  //------------------------------------------------------------------------
  inline
  bool isInside(const RealType* point) const {
     return true;
  }
   
  //------------------------------------------------------------------------
  //! Quantize a floating-point-precision point
  //------------------------------------------------------------------------
  inline
  IntPoint quantize(const RealType* pointIn) const {
    POLY_ASSERT(not mCoordinatesModified);
    return DimensionTraits<Dimension, RealType>::constructPoint(pointIn, &center.x, delta, 0);
  }
  
  //------------------------------------------------------------------------
  //! Dequantize an integer point to floating-point-precision 
  //------------------------------------------------------------------------
  inline
  RealPoint dequantize(const CoordHash* pointIn) const {
    POLY_ASSERT(not mCoordinatesModified);
    RealType p[Dimension];
    for (unsigned j = 0; j < Dimension; ++j) {
      p[j] = center[j] + (RealType(pointIn[j]))*delta;
      // p[j] = center[j] + (RealType(pointIn[j]) + 0.5)*delta;
    }
    return DimensionTraits<Dimension, RealType>::constructPoint(p);
  }

  //------------------------------------------------------------------------
  //! Project a point to the bounding circle
  //------------------------------------------------------------------------
  inline
  RealPoint projectPoint(const RealType* rayPoint,
                         const RealType* rayDirection) const {
    RealPoint result;
    bool test = geometry::rayCircleIntersection(rayPoint, rayDirection, 
                                                &center.x, rinf, 
                                                delta, &result.x);
    POLY_CONTRACT_VAR(test);
    POLY_ASSERT(test);
    return result;
  }

  // // //------------------------------------------------------------------------
  // // //! Check if a point is inside inner bounding box
  // // //------------------------------------------------------------------------
  // // inline
  // // bool isInside(const RealType* point) const {
  // //   bool inside = true;
  // //   for (unsigned j = 0; j != Dimension; ++j)
  // //     inside *= (low[j] <= point[j] and point[j] <= high[j]);
  // //   return inside;
  // // }
   
  // // //------------------------------------------------------------------------
  // // //! Quantize a floating-point-precision point
  // // //! NOTE: If the point is inside the inner box, set the new integer
  // // //!       point's index equal to 1
  // // //------------------------------------------------------------------------
  // // inline
  // // IntPoint quantize(const RealType* pointIn) const {
  // //   POLY_ASSERT(!mCoordinatesModified);
  // //   if (not mOuterBox) {
  // //     return DimensionTraits<Dimension, RealType>::constructPoint(pointIn, &low.x, delta, 0);
  // //   } else {
  // //     if (isInside(pointIn))
  // //       return DimensionTraits<Dimension, RealType>::constructPoint(pointIn, &low.x, delta, 0);
  // //     else
  // //       return DimensionTraits<Dimension, RealType>::constructPoint(pointIn, &olow.x, odelta, 0);
  // //   }
  // // }

  // // //------------------------------------------------------------------------
  // // //! Dequantize an integer point to floating-point-precision 
  // // //------------------------------------------------------------------------
  // // inline
  // // RealPoint dequantize(const CoordHash* pointIn, bool inner=true) const {
  // //   POLY_ASSERT(!mCoordinatesModified);
  // //   RealType p[Dimension];
  // //   RealPoint lo, hi;
  // //   RealType dx;
  // //   if (inner) { lo = low;  hi = high;  dx = delta; } 
  // //   else       { lo = olow; hi = ohigh; dx = odelta;}
  // //   for (unsigned j = 0; j != Dimension; ++j) {
  // //      p[j] = (pointIn[j] == mCoordMax) ? hi[j] : lo[j] + pointIn[j]*dx;
  // //   }
  // //   return DimensionTraits<Dimension, RealType>::constructPoint(p);
  // // }

  // // //------------------------------------------------------------------------
  // // //! Project a point to the bounding circle
  // // //------------------------------------------------------------------------
  // // inline
  // // RealPoint projectPoint(const RealType* rayPoint,
  // //                        const RealType* rayDirection) const {
  // //   RealPoint result;
  // //   RealType radius = (mOuterBox) ? orinf : rinf;
  // //   bool test = geometry::rayCircleIntersection(rayPoint, rayDirection, 
  // //                                               &center.x, radius, 
  // //                                               mDegeneracy, &result.x);
  // //   POLY_ASSERT(test);
  // //   return result;
  // // }


  //------------------------------------------------------------------------
  //! 
  //------------------------------------------------------------------------
  inline
  std::vector<RealPoint> clipToInnerBox(const std::vector<RealPoint> ring) const {
    
  }

};


} //end polytope namespace

#endif
