//----------------------------------------------------------------------------//
// 3D implementation of the convex hull algorithm.
//
// Base on some flash code (!) I found at
// http://en.nicoptere.net/?p=1883
//----------------------------------------------------------------------------//
#ifndef __polytope_convexHull_3d__
#define __polytope_convexHull_3d__

#include <iostream>
#include <vector>
#include <set>

#include "polytope.hh"

namespace polytope {

namespace convexHull_helpers {

template<typename CoordType>
class ConvexHull3d {
public:

  //----------------------------------------------------------------------------
  // An internal 3d point class.
  //----------------------------------------------------------------------------
  class Point {
  public:
    CoordType x, y, z;
    unsigned index;
    Point(const CoordType xi = 0, 
          const CoordType yi = 0, 
          const CoordType zi = 0,
          const unsigned i = std::numeric_limits<unsigned>::max()): x(xi), y(yi), z(zi), index(i) {}
    bool operator<(const Point& rhs) const {
      return (x < rhs.x ? true :
              y < rhs.y ? true :
              z < rhs.z ? true : false);
    }
  };

  //------------------------------------------------------------------------------
  // A fuzzy comparison operator for Point.
  //------------------------------------------------------------------------------
  struct FuzzyPointLessThan {
    CoordType fuzz;
    FuzzyPointLessThan(const CoordType ifuzz = 1): fuzz(ifuzz) {}
    bool operator()(const Point& p1, const Point& p2) {
      return (p2.x - p1.x > fuzz ? true :
              p2.y - p1.y > fuzz ? true : 
              p2.z - p1.z > fuzz ? true : false);
    }
  };

  //----------------------------------------------------------------------------
  // The internal face class.
  //----------------------------------------------------------------------------
  class Face {
  public:
    // Public data.
    static const std::vector<Point>* points;
    int i0, i1, i2;
    double nx, ny, nz;
    
    // Constructor.
    Face(const int ii0, const int ii1, const int ii2): i0(ii0), i1(ii1), i2(ii2) { computePlane(); }
    
    // The plane of the face.
    void computePlane() {
      const Point& p0 = (*points)[i0];
      const Point& p1 = (*points)[i1];
      const Point& p2 = (*points)[i2];
      const double x0 = double(p0.x), y0 = double(p0.y), z0 = double(p0.z);
      const double x1 = double(p1.x), y1 = double(p1.y), z1 = double(p1.z);
      const double x2 = double(p2.x), y2 = double(p2.y), z2 = double(p2.z);
      const double x01 = x1 - x0, y01 = y1 - y0, z01 = z1 - z0;
      const double x02 = x2 - x0, y02 = y2 - y0, z02 = z2 - z0;
      nx = y01*z02 - z01*y02;
      ny = z01*x02 - x01*z02;
      nz = x01*y02 - y01*x02;
      const double norm = sqrt(nx*nx + ny*ny + nz*nz);
      nx /= norm;
      ny /= norm;
      nz /= norm;
    }

    // Test if a point is above the face.
    bool isVisible(const Point& p) const {
      const Point& p0 = (*points)[i0];
      const double dx = double(p.x) - double(p0.x);
      const double dy = double(p.y) - double(p0.y);
      const double dz = double(p.z) - double(p0.z);
      return nx*dx + ny*dy + nz*dz > 0;
    }
	
    // Compute the centroid.
    Point centroid() const {
      const Point& p0 = (*points)[i0];
      const Point& p1 = (*points)[i1];
      const Point& p2 = (*points)[i2];
      const double x0 = double(p0.x), y0 = double(p0.y), z0 = double(p0.z);
      const double x1 = double(p1.x), y1 = double(p1.y), z1 = double(p1.z);
      const double x2 = double(p2.x), y2 = double(p2.y), z2 = double(p2.z);
      return Point( CoordType(( x0 + x1 + x2 ) / 3.0 + 0.5),
                    CoordType(( y0 + y1 + y2 ) / 3.0 + 0.5),
                    CoordType(( z0 + z1 + z2 ) / 3.0 + 0.5));
    }
	
    // Flip the face.
    void flip() {
      std::swap(i0, i1);
      computePlane();
    }

    // Equality.
    bool operator==(const Face& rhs) const {
      return (i0 == rhs.i0) and (i1 == rhs.i1) and (i2 == rhs.i2);
    }

    // Inequality.
    bool operator!=(const Face& rhs) const {
      return not (*this == rhs);
    }
  };
		
  //----------------------------------------------------------------------------
  // Compute the centroid of a face with a point.
  //----------------------------------------------------------------------------
  static Point centroid(const std::vector<Point>& points, const int index, const Face& face) {
    const Point& p = points[index];
    const Point& p0 = points[face.i0];
    const Point& p1 = points[face.i1];
    const Point& p2 = points[face.i2];
    return Point(CoordType(0.25*(double(p.x) + double(p0.x) + double(p1.x) + double(p2.x)) + 0.5),
                 CoordType(0.25*(double(p.y) + double(p0.y) + double(p1.y) + double(p2.y)) + 0.5),
                 CoordType(0.25*(double(p.z) + double(p0.z) + double(p1.z) + double(p2.z)) + 0.5));
  }

  //----------------------------------------------------------------------------
  // The main method -- process the given point cloud and return the hull.
  //----------------------------------------------------------------------------
  static std::vector<std::vector<unsigned> > process(const std::vector<Point>& points) {

    // Pre-conditions.
    ASSERT(points.size() > 3);

    // Initialize the face static info.
    Face::points = &points;

    // Creates a face with the first 3 vertices
    Face face(0, 1, 2);
			
    //this is the center of the tetrahedron, all face should point outwards:
    //they should not be visible to the centroid
    Point v = centroid(points, 3, face);

    if (face.isVisible(v)) face.flip();
			
    Face face0(3, face.i0, face.i1);
    if (face0.isVisible(v)) face0.flip();
			
    Face face1(3, face.i1, face.i2);
    if (face1.isVisible(v)) face1.flip();
			
    Face face2(3, face.i2, face.i0);
    if (face2.isVisible(v)) face2.flip();
			
    //store the tetrahedron faces in the valid faces list
    std::vector<Face> validFaces;
    validFaces.push_back(face);
    validFaces.push_back(face0);
    validFaces.push_back(face1);
    validFaces.push_back(face2);

    //so as we have a convex tetrahedron, we can skip the first 4 points
    for (int i = 4; i < points.size(); ++i) {

      //for each avaiable vertices
      v = points[i];
				
      //checks the point's visibility from all faces
      std::vector<Face> visibleFaces;
      for (typename std::vector<Face>::const_iterator faceItr = validFaces.begin();
           faceItr != validFaces.end();
           ++faceItr) {
        if (faceItr->isVisible(v)) visibleFaces.push_back(*faceItr); 
      }
				
      //the vertex is not visible : it is inside the convex hull, keep on
      if (visibleFaces.size() == 0) continue;
				
      //the vertex is outside the convex hull
      //delete all visible faces from the valid List
      for (typename std::vector<Face>::const_iterator faceItr = visibleFaces.begin();
           faceItr != visibleFaces.end();
           ++faceItr) {
        validFaces.erase(std::remove(validFaces.begin(), validFaces.end(), *faceItr), validFaces.end());
      }        
				
      //special case : only one face is visible
      //it's ok to create 3 faces directly for they won't enclose any other point
      if (visibleFaces.size() == 1) {
        const int i0 = visibleFaces[0].i0;
        const int i1 = visibleFaces[0].i1;
        const int i2 = visibleFaces[0].i2;
        validFaces.push_back(Face(i, i0, i1));
        validFaces.push_back(Face(i, i1, i2));
        validFaces.push_back(Face(i, i2, i0));
        continue;
      }
				
      //creates all possible new faces from the visibleFaces
      std::vector<Face> tmpFaces;
      for (typename std::vector<Face>::const_iterator faceItr = visibleFaces.begin();
           faceItr != visibleFaces.end();
           ++faceItr) {
        tmpFaces.push_back(Face(i, faceItr->i0, faceItr->i1));
        tmpFaces.push_back(Face(i, faceItr->i1, faceItr->i2));
        tmpFaces.push_back(Face(i, faceItr->i2, faceItr->i0));
      }
				
      // Look for faces in the tmp set that have no others in front of them.
      for (typename std::vector<Face>::const_iterator tmpItr = tmpFaces.begin();
           tmpItr != tmpFaces.end();
           ++tmpItr) {
        bool useFace = true;
        typename std::vector<Face>::const_iterator otherItr = tmpFaces.begin();
        while (useFace and otherItr != tmpFaces.end()) {
          if (*otherItr != *tmpItr) {
            useFace = not tmpItr->isVisible(otherItr->centroid());
          }
          ++otherItr;
        }
        if (useFace) validFaces.push_back(*tmpItr);
      }
    }

    // Put together then indices of faces and we're done.
    std::vector<std::vector<unsigned> > result;
    result.reserve(validFaces.size());
    for (typename std::vector<Face>::const_iterator faceItr = validFaces.begin();
         faceItr != validFaces.end();
         ++faceItr) {
      result.push_back(std::vector<unsigned>());
      result.back().push_back(points[faceItr->i0].index);
      result.back().push_back(points[faceItr->i1].index);
      result.back().push_back(points[faceItr->i2].index);
    }
    return result;
  }

};

//------------------------------------------------------------------------------
// Compute the face centroid for a triangular face.
//------------------------------------------------------------------------------
template<typename RealType>
void 
faceCentroid(const std::vector<RealType>& points,
             const std::vector<unsigned>& indices,
             RealType& xc,
             RealType& yc,
             RealType& zc) {
  ASSERT(indices.size() == 3);
  xc = (points[3*indices[0]    ] + points[3*indices[1]    ] + points[3*indices[2]    ])/3.0;
  yc = (points[3*indices[0] + 1] + points[3*indices[1] + 1] + points[3*indices[2] + 1])/3.0;
  zc = (points[3*indices[0] + 2] + points[3*indices[1] + 2] + points[3*indices[2] + 2])/3.0;
}

}

//------------------------------------------------------------------------------
// Helper to output the point class.
//------------------------------------------------------------------------------
template<typename CoordType>
inline
std::ostream&
operator<<(std::ostream& os, const typename convexHull_helpers::ConvexHull3d<CoordType>::Point& value) {
  os << "(" << value.x << " " << value.y << " " << value.z << ")";
  return os;
}

//------------------------------------------------------------------------------
// The 3D convex hull itself.  This is the one users should call -- it forwards
// all work to the worker ConvexHull3d class.
//------------------------------------------------------------------------------
template<typename RealType>
PLC<3, RealType>
convexHull_3d(const std::vector<RealType>& points,
              const RealType* low,
              const RealType& dx) {
  typedef int64_t CoordHash;
  typedef convexHull_helpers::ConvexHull3d<CoordHash>::Point Point;

  // Pre-conditions.
  ASSERT(points.size() % 3 == 0);
  const unsigned n = points.size() / 3;

  unsigned i;

  const RealType& xmin = low[0];
  const RealType& ymin = low[1];
  const RealType& zmin = low[2];

  // Convert the input coordinates to unique integer point types.  Simultaneously we 
  // reduce to the unique set of points.
  typedef std::set<Point, convexHull_helpers::ConvexHull3d<CoordHash>::FuzzyPointLessThan> Set;
  Set pointSet;
  for (i = 0; i != n; ++i) {
    pointSet.insert(Point(CoordHash((points[3*i]     - xmin)/dx + 0.5),
                          CoordHash((points[3*i + 1] - ymin)/dx + 0.5),
                          CoordHash((points[3*i + 2] - zmin)/dx + 0.5), 
                          i));
  }
  ASSERT(pointSet.size() <= n);

  // Extract the unique set of points to a vector.
  std::vector<Point> uniquePoints(pointSet.begin(), pointSet.end());
  ASSERT(uniquePoints.size() == pointSet.size());

  // // Blago!
  // std::cout << "Unique points: " << std::endl;
  // for (unsigned k = 0; k != uniquePoints.size(); ++k) {
  //   std::cout << "   ---> " << uniquePoints[k].x << " "<< uniquePoints[k].y << " " << uniquePoints[k].z << " " << uniquePoints[k].index << std::endl;
  // }
  // // Blago!

  // Dispatch the work 
  const std::vector<std::vector<unsigned> > faces = convexHull_helpers::ConvexHull3d<CoordHash>::process(uniquePoints);

  // Read out the data to the PLC and we're done.
  PLC<3, RealType> plc;
  RealType xc, yc, zc;
  for (i = 0; i != faces.size(); ++i) {
    ASSERT(faces[i].size() == 3);
    plc.facets.push_back(std::vector<int>());
    plc.facets.back().push_back(faces[i][0]);
    plc.facets.back().push_back(faces[i][1]);
    plc.facets.back().push_back(faces[i][2]);
    convexHull_helpers::faceCentroid(points, faces[i], xc, yc, zc);
    // std::cerr << "  -----> " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << " : (" << xc << " " << yc << " " << zc << ")" << std::endl;
  }
  return plc;
}

}

#endif
