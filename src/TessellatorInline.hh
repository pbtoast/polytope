namespace polytope
{

//------------------------------------------------------------------------------
//! This helper method creates a piecewise linear complex (PLC) 
//! representing the bounding box containing the given points and 
//! adds the corners of the bounding box to \a points.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
PLC<Dimension, RealType>
Tessellator<Dimension, RealType>::
boundingBox(std::vector<RealType>& points) const
{
  // Find the minimum and maximum coordinates within the point set and 
  // their point indices.
  RealType Min[Dimension], Max[Dimension];
  for (int d = 0; d < Dimension; ++d)
  {
    Min[d] = FLT_MAX; 
    Max[d] = -FLT_MAX;
  }
  for (size_t i = 0; i < points.size()/Dimension; ++i)
  {
    for (int d = 0; d < Dimension; ++d)
    {
      Min[d] = std::min(Min[d], points[Dimension*i+d]);
      Max[d] = std::max(Max[d], points[Dimension*i+d]);
    }
  }
  // Now create the PLC and add the new generators for the corners.
  return boundingBox(Min, Max, points);
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

}
