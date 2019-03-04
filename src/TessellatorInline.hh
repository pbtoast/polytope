namespace polytope
{

//------------------------------------------------------------------------------
//! This helper method creates a piecewise linear complex (PLC) 
//! representing the bounding box containing the given points and 
//! adds the corners of the bounding box to \a points.
//------------------------------------------------------------------------------
// template<int Dimension, typename RealType>
// inline
// PLC<Dimension, RealType>
// Tessellator<Dimension, RealType>::
// boundingBox(std::vector<RealType>& points) const
// {
//   // Find the minimum and maximum coordinates within the point set and 
//   // their point indices.
//   RealType Min[Dimension], Max[Dimension];
//   for (int d = 0; d < Dimension; ++d)
//   {
//     Min[d] = FLT_MAX; 
//     Max[d] = -FLT_MAX;
//   }
//   for (size_t i = 0; i < points.size()/Dimension; ++i)
//   {
//     for (int d = 0; d < Dimension; ++d)
//     {
//       Min[d] = std::min(Min[d], points[Dimension*i+d]);
//       Max[d] = std::max(Max[d], points[Dimension*i+d]);
//     }
//   }
//   // Now create the PLC and add the new generators for the corners.
//   return boundingBox(Min, Max, points);
// }

//------------------------------------------------------------------------------
//! This helper method returns a normalized collection of points
//! given minimal and maximal values across all dimensions. This maintains
//! aspect ratios during the normalization
//------------------------------------------------------------------------------
template<typename RealType>
inline
std::vector<RealType>
normalizePointsWithSameAspectRatio(const std::vector<RealType>& points,
                                   const RealType xmin,
                                   const RealType xmax)
{
   POLY_ASSERT(xmin < xmax);
   const RealType boxInv = 1.0/std::max(1.0e-30, xmax - xmin);
   const unsigned n = points.size();
   std::vector<RealType> result(n);
   for (unsigned i = 0; i < n; ++i) {
      result[i] = (points[i] - xmin)*boxInv;
      POLY_ASSERT(result[i] >= 0.0 and result[i] <= 1.0);
   }
   return result;
}

//------------------------------------------------------------------------------
//! This helper method scales a collection of normalized points by a common
//! value range to maintain aspect ratios.
//------------------------------------------------------------------------------
template<typename RealType>
inline
void
denormalizePointsWithSameAspectRatio(std::vector<RealType>& normalizedPoints,
                                     const RealType xmin,
                                     const RealType xmax)
{
   POLY_ASSERT(xmin < xmax);
   const RealType L = std::max(1.0e-30, xmax - xmin);
   const unsigned n = normalizedPoints.size();
   for (unsigned i = 0; i < n; ++i) {
      POLY_ASSERT(normalizedPoints[i] >= 0.0 and normalizedPoints[i] <= 1.0);
      normalizedPoints[i] = xmin + L*normalizedPoints[i];
      // POLY_ASSERT2(normalizedPoints[i] >= xmin and normalizedPoints[i] <= xmax,
      //              xmin << " " << normalizedPoints[i] << " " << xmax);
   }
}

//------------------------------------------------------------------------------
//! Return a normalized set of coordinates, also returning the bounding low/high points.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
void
computeCommonBounds(const RealType* points,
                    const int Npoints,
                    const double safetyFactor,
                    RealType& minVal,
                    RealType& maxVal)
{
  minVal =  std::numeric_limits<RealType>::max();
  maxVal = -std::numeric_limits<RealType>::max();       
  RealType low[Dimension], high[Dimension];
  geometry::computeBoundingBox<Dimension, RealType>(points, Npoints, true, low, high);
  for (unsigned j = 0; j != Dimension; ++j) {
    POLY_ASSERT(low[j] <= high[j]);
    minVal = std::min(minVal, low[j] );
    maxVal = std::max(maxVal, high[j]);
    POLY_ASSERT(minVal < maxVal);
  }
  // Add a small safety factor.
  const double dx = safetyFactor*(maxVal - minVal);
  minVal -= dx;
  maxVal += dx;
}


//------------------------------------------------------------------------------
//! Return a normalized set of coordinates, also returning the bounding low/high points.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
std::vector<RealType>
Tessellator<Dimension, RealType>::
computeNormalizedPoints(const std::vector<RealType>& points,
                        const std::vector<RealType>& PLCpoints,
                        const bool computeBounds,
                        RealType* low,
                        RealType* high) const
{
  POLY_ASSERT(points.size() % Dimension == 0);
  if (computeBounds) {
    double low1[Dimension], high1[Dimension];
    geometry::computeBoundingBox<Dimension, RealType>(points, true, low1, high1);
    geometry::computeBoundingBox<Dimension, RealType>(PLCpoints, true, low, high);
    for (unsigned j = 0; j != Dimension; ++j) {
      low[j] = std::min(low[j], low1[j]);
      high[j] = std::max(high[j], high1[j]);
      POLY_ASSERT(low[j] <= high[j]);
    }
    // Add a small safety factor.
    for (unsigned j = 0; j != Dimension; ++j) {
      const double dx = 0.005*(high[j] - low[j]);
      low[j] -= dx;
      high[j] += dx;
    }
  }
  double boxInv[Dimension];
  for (unsigned j = 0; j != Dimension; ++j) boxInv[j] = 1.0/std::max(1e-30, high[j] - low[j]);
  const unsigned n = points.size();
  std::vector<RealType> result(n);
  for (unsigned i = 0; i != n; ++i) {
    const unsigned j = i % Dimension;
    result[i] = (points[i] - low[j])*boxInv[j];
    POLY_ASSERT(result[i] >= 0.0 and result[i] <= 1.0);
  }
  POLY_ASSERT(result.size() == points.size());
  return result;
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// Unbounded case.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
std::vector<unsigned>
Tessellator<Dimension, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     const RealType tol,
                     Tessellation<Dimension, RealType>& mesh) const {

  // Compute the bounding box.
  RealType xmin[Dimension], xmax[Dimension];
  geometry::computeBoundingBox<Dimension, RealType>(points, true, xmin, xmax);

  // Hash the input generators to a unique set.
  std::vector<RealType> uniquePoints;
  std::vector<unsigned> result;
  geometry::uniquePoints<Dimension, RealType>(points, xmin, xmax, tol, uniquePoints, result);

  // Tessellate the unique set, and return the mapping.
  this->tessellate(uniquePoints, mesh);
  return result;
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// Bounded by a box.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
std::vector<unsigned>
Tessellator<Dimension, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     RealType* low,
                     RealType* high,
                     const RealType tol,
                     Tessellation<Dimension, RealType>& mesh) const {

  // Hash the input generators to a unique set.
  std::vector<RealType> uniquePoints;
  std::vector<unsigned> result;
  geometry::uniquePoints<Dimension, RealType>(points, low, high, tol, uniquePoints, result);

  // Tessellate the unique set, and return the mapping.
  this->tessellate(uniquePoints, low, high, mesh);
  return result;
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// PLC case.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
std::vector<unsigned>
Tessellator<Dimension, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     const std::vector<RealType>& PLCpoints,
                     const PLC<Dimension, RealType>& geometry,
                     const RealType tol,
                     Tessellation<Dimension, RealType>& mesh) const {

  // Compute the bounding box.
  RealType xmin[Dimension], xmax[Dimension];
  geometry::computeBoundingBox<Dimension, RealType>(PLCpoints, true, xmin, xmax);

  // Hash the input generators to a unique set.
  std::vector<RealType> uniquePoints;
  std::vector<unsigned> result;
  geometry::uniquePoints<Dimension, RealType>(points, xmin, xmax, tol, uniquePoints, result);

  // Tessellate the unique set, and return the mapping.
  this->tessellate(uniquePoints, PLCpoints, geometry, mesh);
  return result;
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// ReducedPLC case.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
std::vector<unsigned>
Tessellator<Dimension, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     const ReducedPLC<Dimension, RealType>& geometry,
                     const RealType tol,
                     Tessellation<Dimension, RealType>& mesh) const {
  return this->tessellateDegenerate(points, geometry.points, geometry, tol, mesh);
}



//------------------------------------------------------------------------------
// Do a tessellation by first normalizing the input points.
// Unbounded.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
void
Tessellator<Dimension, RealType>::
tessellateNormalized(const std::vector<RealType>& points,
                     Tessellation<Dimension, RealType>& mesh) const {

  // Compute the bounding box.
  RealType xmin, xmax;
  computeCommonBounds<Dimension, RealType>(&points[0], points.size(), 0.0, xmin, xmax);

  // Normalize the input points
  std::vector<RealType> normPoints = 
     normalizePointsWithSameAspectRatio<RealType>(points, xmin, xmax);

  // Tessellate the normalized points
  this->tessellate(normPoints, mesh);

  // Denormalize the mesh nodes
  denormalizePointsWithSameAspectRatio<RealType>(mesh.nodes, xmin, xmax);
}


//------------------------------------------------------------------------------
// Do a tessellation by first normalizing the input points.
// Bounding box.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
void
Tessellator<Dimension, RealType>::
tessellateNormalized(const std::vector<RealType>& points,
                     RealType* low,
                     RealType* high,
                     Tessellation<Dimension, RealType>& mesh) const {

  // Compute the common bounds
  RealType xmin =  std::numeric_limits<RealType>::max();
  RealType xmax = -std::numeric_limits<RealType>::max();      
  for (unsigned j = 0; j < Dimension; ++j) {
    POLY_ASSERT(low[j] <= high[j]);
    xmin = std::min(xmin, low[j] );
    xmax = std::max(xmax, high[j]);
  }
  POLY_ASSERT(xmin < xmax);

  // Normalize the input points
  std::vector<RealType> normPoints = 
     normalizePointsWithSameAspectRatio<RealType>(points, xmin, xmax);
  
  // Normalize the bounds
  RealType normLow[Dimension], normHigh[Dimension];
  for (unsigned j = 0; j < Dimension; ++j) {
    normLow[j]  = (low[j]  - xmin) / (std::max(1.0e-30, xmax - xmin));
    normHigh[j] = (high[j] - xmin) / (std::max(1.0e-30, xmax - xmin));
  }

  // Tessellate the normalized points
  this->tessellate(normPoints, normLow, normHigh, mesh);

  // Denormalize the mesh nodes
  denormalizePointsWithSameAspectRatio<RealType>(mesh.nodes, xmin, xmax);
}

//------------------------------------------------------------------------------
// Do a tessellation by first normalizing the input points.
// PLC.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
void
Tessellator<Dimension, RealType>::
tessellateNormalized(const std::vector<RealType>& points,
                     const std::vector<RealType>& PLCpoints,
                     const PLC<Dimension, RealType>& geometry,
                     Tessellation<Dimension, RealType>& mesh) const {

  // Compute the common bounds
  RealType xmin1, xmax1, xmin2, xmax2;
  computeCommonBounds<Dimension, RealType>(&points[0],    points.size(),    0.0, xmin1, xmax1);
  computeCommonBounds<Dimension, RealType>(&PLCpoints[0], PLCpoints.size(), 0.0, xmin2, xmax2);
  const RealType xmin = std::min(xmin1, xmin2);
  const RealType xmax = std::max(xmax1, xmax2);

  // Normalize the input points
  std::vector<RealType> normPoints = 
     normalizePointsWithSameAspectRatio<RealType>(points, xmin, xmax);

  // Normalize the PLC points
  std::vector<RealType> normPLCPoints = 
     normalizePointsWithSameAspectRatio<RealType>(PLCpoints, xmin, xmax);
  
  // Tessellate the normalized points
  this->tessellate(normPoints, normPLCPoints, geometry, mesh);

  // Denormalize the mesh nodes
  denormalizePointsWithSameAspectRatio<RealType>(mesh.nodes, xmin, xmax);
}

//------------------------------------------------------------------------------
// Do a tessellation by first normalizing the input points.
// ReducedPLC.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
void
Tessellator<Dimension, RealType>::
tessellateNormalized(const std::vector<RealType>& points,
                     const ReducedPLC<Dimension, RealType>& geometry,
                     Tessellation<Dimension, RealType>& mesh) const {
   this->tessellateNormalized(points, geometry.points, geometry, mesh);
}

}
