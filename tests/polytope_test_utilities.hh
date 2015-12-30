//------------------------------------------------------------------------------
// A collection of random stuff useful for testing in polytope.
//------------------------------------------------------------------------------
#ifndef __polytope_test_utilities__
#define __polytope_test_utilities__

#include <sstream>
#include "polytope.hh"
#include "polytope_geometric_utilities.hh"

template<typename RealType> class Boundary2D;

#ifdef HAVE_BOOST
// We use the Boost.Geometry library to handle polygon intersections and such.
// #include <boost/geometry.hpp>
// #include <boost/geometry/geometries/geometries.hpp>
#endif

namespace polytope {

//------------------------------------------------------------------------------
// A macro for checking true/false test conditions.
//------------------------------------------------------------------------------
#define POLY_CHECK(x) { if (!(x)) { std::cout << "FAIL: " << #x << std::endl; exit(-1); } }
#define POLY_CHECK2(x, msg) { if (!(x)) { std::cout << "FAIL: " << #x << std::endl << msg << std::endl; exit(-1); } }

#ifdef HAVE_BOOST
typedef boost::geometry::cs::cartesian cart;
#endif

//------------------------------------------------------------------------------
// Define our own local random number generator wrapping the standard srand &
// rand methods.
//------------------------------------------------------------------------------
inline
double random01() {
  return double(rand())/RAND_MAX;
}


//------------------------------------------------------------------------------
// A simple mesh output function for the SiloWriter
//------------------------------------------------------------------------------
// 2D
template <typename RealType>
void outputMesh(const Tessellation<2,RealType>& mesh,
		std::string prefix,
		const std::vector<RealType>& points,
		const unsigned testCycle = 1,
		const RealType time = 0.0) {
  POLY_ASSERT(points.empty() or
              points.size() == 2*mesh.cells.size());
#ifdef HAVE_SILO
  std::vector<double> index(mesh.cells.size());
  std::vector<double> genx (mesh.cells.size());
  std::vector<double> geny (mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i){
    index[i] = double(i);
    if (!points.empty()) {
      genx[i] = points[2*i  ];
      geny[i] = points[2*i+1];
    }
  }
  std::map<std::string,double*> nodeFields, edgeFields, faceFields, cellFields;
  cellFields["cell_index"] = &index[0];
  cellFields["gen_x"     ] = &genx[0];
  cellFields["gen_y"     ] = &geny[0];
  std::ostringstream os;
  os << prefix;
  polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
					 faceFields, cellFields, os.str(),
					 testCycle, time);
#endif
}

//..............................................................................
// 3D
template <typename RealType>
void outputMesh(const Tessellation<3,RealType>& mesh,
		std::string prefix,
		const std::vector<RealType>& points,
		const unsigned testCycle = 1,
		const RealType time = 0.0) {
  POLY_ASSERT(points.empty() ||
              points.size() == 3*mesh.cells.size());
#ifdef HAVE_SILO
  std::vector<double> index(mesh.cells.size());
  std::vector<double> genx (mesh.cells.size());
  std::vector<double> geny (mesh.cells.size());
  std::vector<double> genz (mesh.cells.size());
  std::vector<double> vol  (mesh.cells.size());
  double cent[3];
  for (int i = 0; i < mesh.cells.size(); ++i){
    index[i] = double(i);
    if (!points.empty()) {
      genx[i] = points[3*i  ];
      geny[i] = points[3*i+1];
      genz[i] = points[3*i+2];
    }
    geometry::computeCellCentroidAndSignedVolume(mesh, i, cent, vol[i]);
  }
  std::map<std::string,double*> nodeFields, edgeFields, faceFields, cellFields;
  cellFields["cell_index"] = &index[0];
  cellFields["gen_x"     ] = &genx[0];
  cellFields["gen_y"     ] = &geny[0];
  cellFields["gen_z"     ] = &genz[0];
  cellFields["volume"    ] = &vol[0];
  std::ostringstream os;
  os << prefix;
  polytope::SiloWriter<3, double>::write(mesh, nodeFields, edgeFields, 
					 faceFields, cellFields, os.str(),
					 testCycle, time);
#endif
}

//------------------------------------------------------------------------------
// Some specialized subsets of outputMesh
//------------------------------------------------------------------------------
template <int nDim, typename RealType>
void outputMesh(const Tessellation<nDim,RealType>& mesh,
		std::string prefix,
		const unsigned testCycle) {
  std::vector<RealType> points;
  outputMesh(mesh, prefix, points, testCycle, 0.0);
}
//------------------------------------------------------------------------------
template <int nDim, typename RealType>
void outputMesh(const Tessellation<nDim,RealType>& mesh,
		std::string prefix) {
  std::vector<RealType> points;
  outputMesh(mesh, prefix, points, 1, 0.0);
}
//------------------------------------------------------------------------------
// a cell-centered field given
template <typename RealType>
void outputMesh(Tessellation<2,RealType>& mesh,
		std::string prefix,
		const std::vector<RealType>& points,
                std::vector<RealType>& cellField,
		const unsigned testCycle = 1,
		const RealType time = 0.0) {
#ifdef HAVE_SILO
  POLY_ASSERT(cellField.size() == mesh.cells.size());
  std::vector<double> index(mesh.cells.size());
  std::vector<double> genx (mesh.cells.size());
  std::vector<double> geny (mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i){
    index[i] = double(i);
    if (!points.empty()) {
      genx[i] = points[2*i  ];
      geny[i] = points[2*i+1];
    }
  }
  std::map<std::string,double*> nodeFields, edgeFields, faceFields, cellFields;
  cellFields["cell_index"] = &index[0];
  cellFields["gen_x"     ] = &genx[0];
  cellFields["gen_y"     ] = &geny[0];
  cellFields["cond"      ] = &cellField[0];
  std::ostringstream os;
  os << prefix;
  polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
					 faceFields, cellFields, os.str(),
					 testCycle, time);
#endif
}



//------------------------------------------------------------------------------
// Wrapper to tessellate a 2D boundary for both VoroPP_2D and Triangle
//------------------------------------------------------------------------------
template <typename RealType>
void tessellate2D(std::vector<RealType>& points,
                  Boundary2D<RealType>& boundary,
                  Tessellator<2,RealType>& tessellator,
                  Tessellation<2,RealType>& mesh){
   if( tessellator.handlesPLCs() ){
      tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
   }else{
      tessellator.tessellate(points, boundary.mLow, boundary.mHigh, mesh);
   }
}

#ifdef HAVE_BOOST
//------------------------------------------------------------------------------
// Make a 2D Boost.Geometry point
//------------------------------------------------------------------------------
template <typename RealType>
boost::geometry::model::point<RealType,2,cart> makePoint2D(std::vector<RealType>& position) {
   POLY_ASSERT( position.size() == 2 );
   boost::geometry::model::point<RealType,2,cart> point;
   boost::geometry::set<0>(point, position[0]);
   boost::geometry::set<1>(point, position[1]);
   return point;
}

//------------------------------------------------------------------------------
// Make a 3D Boost.Geometry point
//------------------------------------------------------------------------------
template <typename RealType>
boost::geometry::model::point<RealType,3,cart> makePoint3D(std::vector<RealType>& position) {
   POLY_ASSERT( position.size() == 3 );
   boost::geometry::model::point<RealType,3,cart> point;
   boost::geometry::set<0>(point, position[0]);
   boost::geometry::set<1>(point, position[1]);
   boost::geometry::set<2>(point, position[2]);
   return point;
}

//------------------------------------------------------------------------------
// Make a Boost.Geometry polygon from a concatenated vector of (x,y) points
//------------------------------------------------------------------------------
template <typename RealType>
boost::geometry::model::polygon<boost::geometry::model::point<RealType,2,cart>,false> 
makePolygon( std::vector<RealType>& points ) {
   typedef boost::geometry::model::point<RealType,2,cart> RealPoint;
   boost::geometry::model::polygon<RealPoint,false> polygon;
   for (unsigned i = 0; i < points.size()/2; ++i) {
      boost::geometry::append( polygon, RealPoint(points[2*i],points[2*i+1]) );
   }
   boost::geometry::append( polygon, RealPoint(points[0],points[1]) );
   return polygon;
}

//------------------------------------------------------------------------------
// Make a Boost.Geometry polygon from a PLC and its point list
//------------------------------------------------------------------------------
template <typename RealType>
boost::geometry::model::polygon<boost::geometry::model::point<RealType,2,cart>,false> 
makePolygon( PLC<2,RealType>& PLC, std::vector<RealType>& PLCpoints ) {
   typedef boost::geometry::model::point<RealType,2,cart> RealPoint;
   typedef boost::geometry::model::polygon<RealPoint,false> BGpolygon;
   
   unsigned i,j;
   BGpolygon polygon;
   // Walk the facets and add the first node
   for (j = 0; j != PLC.facets.size(); ++j){
      POLY_ASSERT( PLC.facets[j].size() == 2 );
      i = PLC.facets[j][0];
      boost::geometry::append( polygon, RealPoint(PLCpoints[2*i],PLCpoints[2*i+1]) );
   }
   i = PLC.facets[0][0];
   boost::geometry::append( polygon, RealPoint(PLCpoints[2*i],PLCpoints[2*i+1]) ); //Close the polygon
   
   // Walk the facets composing each hole and add the first node
   const unsigned nHoles = PLC.holes.size();
   if (nHoles > 0) {
      typename BGpolygon::inner_container_type& holes = polygon.inners();
      holes.resize(nHoles);
      for (unsigned ihole = 0; ihole != nHoles; ++ihole){
         for (j = 0; j != PLC.holes[ihole].size(); ++j){
            POLY_ASSERT( PLC.holes[ihole][j].size() == 2 );
            i = PLC.holes[ihole][j][0];
            boost::geometry::append( holes[ihole], RealPoint( PLCpoints[2*i], PLCpoints[2*i+1] ) );
         }
         i = PLC.holes[ihole][0][0];
         boost::geometry::append( holes[ihole], RealPoint( PLCpoints[2*i], PLCpoints[2*i+1] ) );  //Close the polygon
      }
   }
   return polygon;
}

//------------------------------------------------------------------------------
// Compute the area of a polytope tessellation cell-by-cell using Boost.Geometry
//------------------------------------------------------------------------------
template <typename RealType>
RealType computeTessellationArea( polytope::Tessellation<2,RealType>& mesh ) {
   typedef boost::geometry::model::point<RealType,2,cart> RealPoint;
   RealType area = 0;
   for (unsigned i = 0; i != mesh.cells.size(); ++i) {
      std::vector<RealType> nodeCell;
      for (std::vector<int>::const_iterator faceItr = mesh.cells[i].begin();
           faceItr != mesh.cells[i].end(); ++faceItr){
         const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
         POLY_ASSERT(iface < mesh.faceCells.size());
         POLY_ASSERT(mesh.faces[iface].size() == 2);
         const unsigned inode = *faceItr < 0 ? mesh.faces[iface][1] : mesh.faces[iface][0];
         nodeCell.push_back( mesh.nodes[2*inode  ] );
         nodeCell.push_back( mesh.nodes[2*inode+1] );
      }
      boost::geometry::model::polygon<RealPoint,false> cellPolygon = makePolygon( nodeCell );
      area += boost::geometry::area( cellPolygon );
      nodeCell.clear();
   }
   return area;
}
#endif

}

#endif
