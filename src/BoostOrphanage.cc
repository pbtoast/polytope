//------------------------------------------------------------------------
// BoostOrphanage
//------------------------------------------------------------------------
#include <set>
#include <map>
#include <vector>

#include "polytope.hh"

namespace polytope {

using namespace std;


//------------------------------------------------------------------------
// A collection of helper functions
//------------------------------------------------------------------------

namespace {

//------------------------------------------------------------------------
// Union a Boost.Geometry ring with a Boost.Geometry multi_polygon.
// The resulting multi_polygon is corrected to ensure it conforms to
// the proper geometric concept.
//------------------------------------------------------------------------
void
createBGUnion(boost::geometry::model::ring<Point2<int64_t>,false> ring, 
              boost::geometry::model::multi_polygon<
                boost::geometry::model::polygon<Point2<int64_t>,false> >& multiPolygon ) {
  boost::geometry::model::multi_polygon<
    boost::geometry::model::polygon<Point2<int64_t>,false> > temp;
  boost::geometry::union_(multiPolygon, ring, temp);
  boost::geometry::correct(temp);
  multiPolygon=temp;
}

} // end anonymous namespace




//------------------------------------------------------------------------------
template<typename RealType>
BoostOrphanage<RealType>::
BoostOrphanage(const Tessellator<2, RealType>* tessellator):
  OrphanageBase<2, RealType>(tessellator) {
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
BoostOrphanage<RealType>::
~BoostOrphanage() {
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
BoostOrphanage<RealType>::
adoptOrphans(const vector<RealType>& points,
	     const RealType* low,
	     const RealType* high,
	     const RealType dx,
	     vector<BGring>& cellRings,
	     vector<BGring>& orphans) const {
  // Pre-conditions
  POLY_ASSERT(!points.empty());
  POLY_ASSERT(!cellRings.empty());
  POLY_ASSERT(points.size()/2 == cellRings.size());
  POLY_ASSERT(!orphans.empty());
  
  const int numGenerators = points.size()/2;
  int i;
  
  // Construct map from node points to neighboring cells
  std::map<IntPoint, std::set<int> > point2neighbors;
  for (i = 0; i != cellRings.size(); ++i){
    for (typename BGring::const_iterator itr = cellRings[i].begin();
         itr != cellRings[i].end() - 1; ++itr){
      point2neighbors[*itr].insert(i);
    }
  }
  
  // Create reduced set of orphans by unioning all that neighbor each other.
  // A Boost.Geometry multi-polygon does this operation seemlessly
  BGmulti_polygon orphanUnion;
  for (i = 0; i != orphans.size(); ++i) {
    createBGUnion(orphans[i], orphanUnion);
  }
  
  // Compute the adoption on each orphan in the reduced set
  for (i = 0; i < orphanUnion.size(); ++i) {
    BGring orphan = orphanUnion[i].outer();
    POLY_ASSERT(!orphan.empty());
    POLY_ASSERT(orphan.front() == orphan.back());
    
    // Determine neighbors
    std::set<int> orphanNeighbors;
    for (typename BGring::const_iterator pointItr = orphan.begin();
         pointItr != orphan.end()-1; ++pointItr) {
      std::map<IntPoint, std::set<int> >::iterator it = point2neighbors.find(*pointItr);
      if (it != point2neighbors.end()){
        std::set<int> neighborSet = it->second;
        orphanNeighbors.insert(neighborSet.begin(), neighborSet.end());
      }
    }
    POLY_ASSERT( orphanNeighbors.size() > 0 );

    // If the orphan only has a single neighbor, we can skip a lot of work.
    // No need to tessellate - simply union the orphan with its neighbor cell.
    std::vector<BGring> subCellRings;
    if (orphanNeighbors.size() > 1) {
       
      // Union the orphan and its neighbors to get a bounded polygon
      BGmulti_polygon neighborCells;
      createBGUnion(orphan,neighborCells);
      
      // Collect the neighboring generators
      std::vector<RealType> subpoints;
      for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
           nbItr != orphanNeighbors.end(); ++nbItr){
        subpoints.push_back( points[2*(*nbItr)  ] );
        subpoints.push_back( points[2*(*nbItr)+1] );
        createBGUnion(cellRings[*nbItr],neighborCells);
      }
      POLY_ASSERT2( neighborCells.size() > 0, "Union produced empty set!" );
      
      // Throw an error if the bounding polygon has more than one piece
      if (neighborCells.size() > 1) {
        std::cerr << "Blago!" << std::endl;
        for (i = 0; i != neighborCells.size(); ++i){
          std::cerr << "Polygon " << i << " in the union contains" << std::endl;
          for (typename BGring::const_iterator itr = neighborCells[i].outer().begin();
               itr != neighborCells[i].outer().end(); ++itr){
            std::cerr << (*itr)
                      << "(" << (*itr).realx(low[0], dx) 
                      << "," << (*itr).realy(low[1], dx) << ")" << std::endl;
          }
          POLY_ASSERT(false);
        }
      }
      BGring boundaryRing = neighborCells[0].outer();
      POLY_ASSERT(boundaryRing.front() == boundaryRing.back());
      POLY_ASSERT(boundaryRing.size() > 2);
      
      // TODO: Make sure union-ing rings that share a common face results in 
      //       a closed boundary, has no repeated nodes, etc. etc.
      
      // TODO: Check whether converting the PLC points back to doubles to compute
      //       the sub-tessellation gives a valid full tessellation after the
      //       cell adoption loop concludes
      
      // Extract the vertices of the bounding polygon
      std::vector<RealType> subPLCpoints;
      int nSides = 0;
      for (typename BGring::const_iterator itr = boundaryRing.begin();
           itr != boundaryRing.end() - 1; ++itr, ++nSides) {
        subPLCpoints.push_back( (*itr).realx(low[0], dx) );
        subPLCpoints.push_back( (*itr).realy(low[1], dx) );
      }
      
      // Form the bounding PLC
      PLC<2, RealType> subPLC;
      subPLC.facets.resize(nSides, std::vector<int>(2) );
      for (unsigned ii = 0; ii < nSides; ++ii) {
        subPLC.facets[ii][0] = ii;
        subPLC.facets[ii][1] = (ii+1) % nSides;
      }
      
      // Tessellate this sub-region
      Tessellation<2,RealType> submesh;
      this->callPrivateTessellate(subpoints, subPLCpoints, subPLC, 
                                  low, high, dx, submesh);
      
      // Convert the tessellation cells to integer cell rings
      convertTessellationToRings(submesh, &low[0], dx, subCellRings);
      
      // // Blago!
      // std::cerr << std::endl << "Subpoints:" << std::endl;
      // for (int ii = 0; ii != subpoints.size(); ++ii){
      //   std::cerr << "  (" << subpoints[2*ii] << "," 
      //             << subpoints[2*ii+1] << ")" << std::endl;
      // }
      // std::cerr << std::endl << "SubPLCpoints:" << std::endl;
      // for (int ii = 0; ii != subPLCpoints.size(); ++ii){
      //   std::cerr << "  (" << subPLCpoints[2*ii] << "," 
      //             << subPLCpoints[2*ii+1] << ")" << std::endl;
      // }
      // std::cerr << std::endl << "SubCellRings:" << std::endl;
      // for (int ii = 0; ii != subCellRings.size(); ++ii) {
      //   std::cerr << "  Cell " << ii << std::endl;
      //   for (typename BGring::const_iterator itr = subCellRings[i].begin();
      //        itr != subCellRings[i].end(); ++itr){
      //     std::cerr << (*itr).realx(low[0], dx) << " " 
      //               << (*itr).realy(low[1], dx) << " ";
      //   }
      // }
      // // Blago!
      
    }
    
    // We're only concerned with the cells in the sub-tessellation whose generators
    // are immediate neighbors of the orphaned chunk. These are the only cells which can
    // "adopt" the orphan based on a local Voronoi principle
    for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
         nbItr != orphanNeighbors.end(); ++nbItr){
      std::set<int>::iterator it = orphanNeighbors.find(*nbItr);
      POLY_ASSERT(it != orphanNeighbors.end());
      // Index for the sub-tessellation
      int subIndex = distance(orphanNeighbors.begin(), it);  
      POLY_ASSERT(subIndex < orphanNeighbors.size());
      // Index for the full tessellation
      int thisIndex = *nbItr;                                
      POLY_ASSERT(thisIndex < numGenerators);
      // The full tessellation's cell ring
      BGring thisRing;
      
      if (orphanNeighbors.size() > 1) {
        thisRing = subCellRings[subIndex];
        
        // Simplify the resulting ring. Removes points that are within some minimum
        // distance to their neighbors. Setting distance = 2 merges ring elements
        // that are within one quantized mesh spacing. This essentially removes
        // repeated cell nodes having length-zero cell faces.
        BGring simplifiedRing;
        boost::geometry::simplify(thisRing, simplifiedRing, 2);
        thisRing = simplifiedRing;
        POLY_ASSERT(thisRing.size() > 2);
        POLY_ASSERT(thisRing.front() == thisRing.back());
      }	
      
      // If the orphan has only a single neighbor, just compute its union with
      // that neighbor's cell ring from the full tessellation
      else{
        thisRing = orphan;
      }
      
      // Union this new cell ring with the original cell ring from the full tessellation
      std::vector<BGring> unionRing;
      boost::geometry::union_(thisRing, cellRings[thisIndex], unionRing);
      if (unionRing.size() > 1) {
        std::cerr << "Blago!" << std::endl << "Cell " << thisIndex
                  << " has more than one cell ring:" << std::endl;
        for( i=0; i<unionRing.size(); ++i){
          std::cerr << std::endl << "Ring " << i << ":" << std::endl;
          for (typename BGring::const_iterator itr = unionRing[i].begin();
               itr != unionRing[i].end(); ++itr){
            std::cerr << (*itr).realx(low[0], dx) << " " 
                      << (*itr).realy(low[1], dx) << " "
                      << (*itr) << std::endl;
          }
        }
      }
      POLY_ASSERT(unionRing.size() == 1);
      thisRing = unionRing[0];
      
      // Simplify the final ring. 
      BGring simplifiedRing;
      boost::geometry::simplify(thisRing, simplifiedRing, 1);
      thisRing = simplifiedRing;
      POLY_ASSERT(thisRing.size() > 2);
      POLY_ASSERT(thisRing.front() == thisRing.back());
      
      // // Blago!
      // std::cerr << "\nCell " << thisIndex << " before final cleaning" << std::endl;
      // for (typename BGring::const_iterator itr = thisRing.begin();
      //      itr != thisRing.end(); ++itr){
      //    std::cerr << (*itr).realx(low[0], dx) << " " 
      //              << (*itr).realy(low[1], dx) << " "
      //              << (*itr) << std::endl;
      // }
      // // Blago!

      // Union may produce cell rings containing faces broken into multiple
      // pieces by redundant hanging nodes. This manifests as three (or more)
      // collinear nodes. We clip them as a final connection step
      BGring finalRing;
      bool collinear;
      for (typename BGring::const_iterator itr = thisRing.begin()+1;
           itr != thisRing.end()-1; ++itr) {
        collinear = geometry::collinear<2,CoordHash>(&(*(itr-1)).x,
                                                     &(*(itr  )).x,
                                                     &(*(itr+1)).x,
                                                     1);
        if(!collinear) boost::geometry::append(finalRing,*itr);
      }
      
      // Check the beginning/end point
      collinear = geometry::collinear<2,CoordHash>(&(*(thisRing.end()-2  )).x,
                                                   &(*(thisRing.begin()  )).x,
                                                   &(*(thisRing.begin()+1)).x,
                                                   1);
      if(!collinear) boost::geometry::append(finalRing,*(thisRing.begin()));

      // Close the ring
      boost::geometry::correct(finalRing);
      POLY_ASSERT(finalRing.size() > 2);
      POLY_ASSERT(finalRing.front() == finalRing.back());
      thisRing = finalRing;


      // // Blago!
      // std::cerr << "\nCell " << thisIndex << " after final cleaning" << std::endl;
      // for (typename BGring::const_iterator itr = thisRing.begin();
      //      itr != thisRing.end(); ++itr){
      //    std::cerr << (*itr).realx(low[0], dx) << " " 
      //              << (*itr).realy(low[1], dx) << " "
      //              << (*itr) << std::endl;
      // }
      // // Blago!
      
      // Replace the corresponding cell ring
      cellRings[thisIndex] = thisRing;
    }
  }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Explicit instantiation
//------------------------------------------------------------------------------
template class BoostOrphanage<double>;

} //end polytope namespace
