//------------------------------------------------------------------------
// BoostOrphanage
//------------------------------------------------------------------------
#include <set>
#include <map>
#include <vector>

#include "polytope.hh"
#include "BoostOrphanage.hh"

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

//------------------------------------------------------------------------
// Union may produce cell rings containing faces broken into multiple
// pieces by redundant hanging nodes. This manifests as three (or more)
// collinear nodes. We clip them as a final connection step
//------------------------------------------------------------------------
void
cleanRingEdges(boost::geometry::model::ring<Point2<int64_t>,false>& ring) {
  // Pre-conditions
  POLY_ASSERT(ring.size() > 2);
  POLY_ASSERT(ring.front() == ring.back());
  //POLY_ASSERT(!boost::geometry::intersects(ring));

  // Initialize the temporary ring
  typedef boost::geometry::model::ring<Point2<int64_t>,false> BGring;
  BGring tmpRing;
  
  // Check collinearity on the middle points
  bool collinear;
  for (BGring::const_iterator itr = ring.begin()+1;
       itr != ring.end()-1; ++itr) {
    collinear = geometry::collinear<2,CoordHash>(&(*(itr-1)).x,
                                                 &(*(itr  )).x,
                                                 &(*(itr+1)).x,
                                                 1);
    if(!collinear) boost::geometry::append(tmpRing,*itr);
  }
  
  // Check the beginning/end point
  collinear = geometry::collinear<2,CoordHash>(&(*(ring.end()-2  )).x,
                                               &(*(ring.begin()  )).x,
                                               &(*(ring.begin()+1)).x,
                                               1);
  if(!collinear) boost::geometry::append(tmpRing,*(ring.begin()));

  // Close the ring
  boost::geometry::correct(tmpRing);

  // Post-conditions
  POLY_ASSERT(tmpRing.size() > 2);
  POLY_ASSERT(tmpRing.front() == tmpRing.back());
  POLY_ASSERT(!boost::geometry::intersects(tmpRing));
  ring = tmpRing;
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
	     const QuantizedCoordinates<2, RealType>& coords,
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
    POLY_ASSERT2(!boost::geometry::intersects(orphan),
                 "Unioning orphans has produced self-intersections");
    
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
               itr != neighborCells[i].outer().end(); 
               ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
          POLY_ASSERT(false);
        }
      }
      BGring boundaryRing = neighborCells[0].outer();
      BGpolygon boundaryPolygon = neighborCells[0];
      {
        BGring simplifiedRing;
        boost::geometry::simplify(boundaryRing, simplifiedRing, 3);
        boundaryRing = simplifiedRing;
        cleanRingEdges(boundaryRing);
      }
      
      // // Blago!
      // std::cerr << "\nBoundary ring" << std::endl;
      // for (typename BGring::const_iterator itr = boundaryRing.begin();
      //      itr != boundaryRing.end();
      //      ++itr)   std::cerr << (*itr) << coords.dequantize(*itr) << std::endl;
      // std::cerr << std::endl;
      // // Blago!


      // TODO: Check whether converting the PLC points back to doubles to compute
      //       the sub-tessellation gives a valid full tessellation after the
      //       cell adoption loop concludes
      // SPOILER ALERT! It doesn't!


      // Extract the vertices of the bounding polygon
      std::vector<CoordHash> subIntPLCpoints;
      int nSides = 0;
      for (typename BGring::const_iterator itr = boundaryRing.begin();
           itr != boundaryRing.end() - 1; ++itr, ++nSides) {
         subIntPLCpoints.push_back((*itr).x);
         subIntPLCpoints.push_back((*itr).y);
      }
      
      // Form the bounding PLC
      PLC<2, RealType> subPLC;
      subPLC.facets.resize(nSides, std::vector<int>(2) );
      for (unsigned ii = 0; ii < nSides; ++ii) {
        subPLC.facets[ii][0] = ii;
        subPLC.facets[ii][1] = (ii+1) % nSides;
      }
      
      // Tessellate this sub-region
      // Tessellation<2,RealType> submesh;
      // this->callPrivateTessellate(subpoints, subIntPLCpoints, subPLC, coords, submesh);
      vector<vector<vector<CoordHash> > > IntCells;
      this->callPrivateTessellate(subpoints, subIntPLCpoints, subPLC, coords, IntCells);
      
      subCellRings.resize(IntCells.size());
      for (unsigned ii = 0; ii != IntCells.size(); ++ii) {
        vector<IntPoint> cellNodes;
        for (unsigned jj = 0; jj != IntCells[ii].size(); ++jj) {
          POLY_ASSERT(IntCells[ii][jj].size() == 2);
          cellNodes.push_back(IntPoint(IntCells[ii][jj][0], IntCells[ii][jj][1]));
        }
        boost::geometry::assign(subCellRings[ii], BGring(cellNodes.begin(),
                                                         cellNodes.end()));
        POLY_ASSERT(IntCells[ii].size()      == subCellRings[ii].size());
        POLY_ASSERT(subCellRings[ii].front() == subCellRings[ii].back());
      }


      // // Blago!
      // {
      //   vector<double> px(submesh.cells.size());
      //   vector<double> py(submesh.cells.size());
      //   for (unsigned i = 0; i != submesh.cells.size(); ++i) {
      //     px[i] = subpoints[2*i];
      //     py[i] = subpoints[2*i+1];
      //   }
      //   map<string, double*> fields, cellFields;
      //   cellFields["gen_x"] = &px[0];
      //   cellFields["gen_y"] = &py[0];
      //   SiloWriter<2, RealType>::write(submesh, fields, fields, 
      //                                  fields, cellFields, "test_BoostOrphanage");
      // }
      // std::cerr << std::endl << "Subpoints:" << std::endl;
      // for (int ii = 0; ii != subpoints.size()/2; ++ii){
      //   std::cerr << "  (" << subpoints[2*ii] << "," 
      //             << subpoints[2*ii+1] << ")" << std::endl;
      // }
      // std::cerr << std::endl << "SubPLC:" << std::endl;
      // std::cerr << subPLC << std::endl;
      // std::cerr << std::endl << "SubPLCpoints:" << std::endl;
      // for (int ii = 0; ii != subIntPLCpoints.size()/2; ++ii){
      //    std::cerr << coords.dequantize(IntPoint(subIntPLCpoints[2*ii],subIntPLCpoints[2*ii+1])) << std::endl;
      // }
      // // Blago!


      // Convert the tessellation cells to integer cell rings
      // convertTessellationToRings(submesh, &coords.low[0], coords.delta, subCellRings);



      // // Blago!
      // std::cerr << std::endl << "SubCellRings:" << std::endl;
      // for (int ii = 0; ii != subCellRings.size(); ++ii) {
      //   std::cerr << "  Cell " << ii << std::endl;
      //   for (typename BGring::const_iterator itr = subCellRings[ii].begin();
      //        itr != subCellRings[ii].end(); 
      //        ++itr)   std::cerr << (*itr) << coords.dequantize(*itr) << std::endl;
      //   std::cerr << std::endl;
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
        int nsimp = 2;
        while(boost::geometry::intersects(thisRing) and nsimp < 14) {
          BGring simplifiedRing;
          boost::geometry::simplify(thisRing, simplifiedRing, nsimp);
          thisRing = simplifiedRing;
          ++nsimp;
        }
        POLY_ASSERT(thisRing.size() > 2);
        POLY_ASSERT(thisRing.front() == thisRing.back());
        POLY_ASSERT2(!boost::geometry::intersects(thisRing),
                     "Subcell " << subIndex << " has self-intersections");
      }
      
      // If the orphan has only a single neighbor, just compute its union with
      // that neighbor's cell ring from the full tessellation
      else{
        thisRing = orphan;
      }
      
      // // Blago!
      // std::cerr << "\nCell " << thisIndex << " before union" << std::endl;
      // for (typename BGring::const_iterator itr = thisRing.begin();
      //      itr != thisRing.end(); 
      //      ++itr)   std::cerr << (*itr) << coords.dequantize(*itr) << std::endl;
      // // Blago!

      // Union this new cell ring with the original cell ring from the full tessellation
      std::vector<BGring> unionRing;
      boost::geometry::union_(thisRing, cellRings[thisIndex], unionRing);
      if (unionRing.size() > 1) {
        std::cerr << "Blago!" << std::endl << "Cell " << thisIndex
                  << " has more than one cell ring:" << std::endl;
        for( i=0; i<unionRing.size(); ++i){
          std::cerr << std::endl << "Ring " << i << ":" << std::endl;
          for (typename BGring::const_iterator itr = unionRing[i].begin();
               itr != unionRing[i].end(); 
               ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
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
      //      itr != thisRing.end(); 
      //      ++itr)   std::cerr << (*itr) << coords.dequantize(*itr) << std::endl;
      // // Blago!

      // Union may produce cell rings containing faces broken into multiple
      // pieces by redundant hanging nodes. This manifests as three (or more)
      // collinear nodes. We clip them as a final connection step
      cleanRingEdges(thisRing);

      // // Blago!
      // std::cerr << "\nCell " << thisIndex << " after final cleaning" << std::endl;
      // for (typename BGring::const_iterator itr = thisRing.begin();
      //      itr != thisRing.end();
      //      ++itr)   std::cerr << (*itr) << coords.dequantize(*itr) << std::endl;
      // // Blago!
      
      // Replace the corresponding cell ring
      cellRings[thisIndex] = thisRing;
    }
  }
}
//------------------------------------------------------------------------------


// //------------------------------------------------------------------------------
// template<typename RealType>
// void
// BoostOrphanage<RealType>::
// callPrivateTessellate(const std::vector<RealType>& points,
//                       const BGpolygon boundary,
//                       const QuantizedCoordinates<2, RealType>& coords,
//                       std::vector<BGring>& cellRings) const {
//   POLY_ASSERT2(false, "Doesn't work!");
// }


//------------------------------------------------------------------------------
// Explicit instantiation
//------------------------------------------------------------------------------
template class BoostOrphanage<double>;

} //end polytope namespace
