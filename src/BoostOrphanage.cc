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
  typedef int64_t CoordHash;
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
  //POLY_ASSERT(!boost::geometry::intersects(tmpRing));
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
  int i,j;
  
  // Create reduced set of orphans by unioning all that neighbor each other.
  // A Boost.Geometry multi-polygon does this operation seemlessly
  BGmulti_polygon reducedOrphanSet;
  for (i = 0; i != orphans.size(); ++i) {
    createBGUnion(orphans[i], reducedOrphanSet);
  }  
  
  // std::cerr << "Number of orphans = " << reducedOrphanSet.size() << std::endl;

  // Compute the adoption on each orphan in the reduced set
  for (i = 0; i < reducedOrphanSet.size(); ++i) {
    BGring orphan = reducedOrphanSet[i].outer();
    POLY_ASSERT(!orphan.empty());
    POLY_ASSERT(orphan.front() == orphan.back());
    POLY_ASSERT2(!boost::geometry::intersects(orphan),
                 "Unioning orphans has produced self-intersections");
    
    // // Blago!
    // std::cerr << "Orphan " << i << ":" << std::endl;
    // for (typename BGring::const_iterator itr = orphan.begin();
    //      itr != orphan.end(); 
    //      ++itr)  std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
    // std::cerr << std::endl;
    // // Blago!

    // Construct map from node points to neighboring cells
    std::map<IntPoint, std::set<int> > point2neighbors;
    for (j = 0; j != cellRings.size(); ++j) {
      for (typename BGring::const_iterator itr = cellRings[j].begin();
	   itr != cellRings[j].end() - 1; ++itr) {
	point2neighbors[*itr].insert(j);
      }
    }
    
    // Determine neighbors
    std::set<int> orphanNeighbors;
    for (typename BGring::const_iterator pointItr = orphan.begin();
         pointItr != orphan.end()-1; 
         ++pointItr) {
      std::map<IntPoint, std::set<int> >::iterator it = 
         point2neighbors.find(*pointItr);
      if (it != point2neighbors.end()) {
        std::set<int> neighborSet = it->second;
        orphanNeighbors.insert(neighborSet.begin(), neighborSet.end());
      }
    }
    POLY_ASSERT(orphanNeighbors.size() > 0);
    
    // If the orphan only has a single neighbor, we can skip a lot of work.
    // No need to tessellate - simply union the orphan with its neighbor cell.
    std::vector<BGring> subCellRings;
    if (orphanNeighbors.size() > 1) {
      
      // Union the orphan and its neighbors to get a bounded polygon
      BGmulti_polygon neighborCells;
      createBGUnion(orphan, neighborCells);
      POLY_ASSERT(neighborCells.size() == 1);
      
      // // Blago!
      // std::cerr << "Initial union polygon" << std::endl;
      // for (typename BGring::const_iterator itr = neighborCells[0].outer().begin();
      //      itr != neighborCells[0].outer().end(); 
      //      ++itr)  std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
      // std::cerr << std::endl;
      // // Blago!

      // Collect the neighboring generators
      std::vector<RealType> subpoints;
      for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
           nbItr != orphanNeighbors.end(); 
           ++nbItr){
        subpoints.push_back( points[2*(*nbItr)  ] );
        subpoints.push_back( points[2*(*nbItr)+1] );
        createBGUnion(cellRings[*nbItr],neighborCells);
      }
      POLY_ASSERT2( neighborCells.size() > 0, "Union produced empty set!" );

      // std::cerr << "Number of neighbors = " << orphanNeighbors.size() << std::endl;
      
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

      // // Blago!
      // std::cerr << "Bounding Union Polygon" << std::endl;
      // for (typename BGring::const_iterator itr = neighborCells[0].outer().begin();
      //      itr != neighborCells[0].outer().end(); 
      //      ++itr)  std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
      // std::cerr << std::endl;
      // // Blago!
      
#ifndef NDEBUG
      for (typename BGring::const_iterator itr = neighborCells[0].outer().begin();
	   itr != neighborCells[0].outer().end() - 1; 
	   ++itr ) {
	bool result = false;
	for (typename BGring::const_iterator pItr = orphan.begin();
	     pItr != orphan.end() - 1;
	     ++pItr) result += (*itr == *pItr);
	for (std::set<int>::const_iterator nbItr = orphanNeighbors.begin();
	     nbItr != orphanNeighbors.end(); 
	     ++nbItr) {
	  for (typename BGring::const_iterator pItr = cellRings[*nbItr].begin();
	       pItr != cellRings[*nbItr].end() - 1;
	       ++pItr) result += (*itr == *pItr);
	}
	POLY_ASSERT2(result, "Union error: there are points in the bounding ring "
		     << "that are not in the orphan or any of its neighbors.");
      }
#endif
      
      BGring boundaryRing = neighborCells[0].outer();
      //BGpolygon boundaryPolygon = neighborCells[0];
      cleanRingEdges(boundaryRing);
      POLY_ASSERT(boundaryRing.size()  >  2);
      POLY_ASSERT(boundaryRing.front() == boundaryRing.back());
      POLY_ASSERT(!boost::geometry::intersects(boundaryRing));

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
      // std::cerr << std::endl << "SubCellRings:" << std::endl;
      // for (int ii = 0; ii != subCellRings.size(); ++ii) {
      //   std::cerr << "  Cell " << ii << std::endl;
      //   for (typename BGring::const_iterator itr = subCellRings[ii].begin();
      //        itr != subCellRings[ii].end(); 
      //        ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
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
        
        // // Blago!
        // std::cerr << "\nCell " << thisIndex << " before fixup" << std::endl;
        // for (typename BGring::const_iterator itr = thisRing.begin();
        //      itr != thisRing.end();
        //      ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
        // std::cerr << "\nCurrent Cell " << thisIndex << std::endl;
        // for (typename BGring::const_iterator itr = cellRings[thisIndex].begin();
        //      itr != cellRings[thisIndex].end();
        //      ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
        // // Blago!


        // When calling the private tessellate routine, we may end up with new cell
        // rings that vertices that are close to the original ones, but off by a few
        // quantized grid spacings. Unioning the old and new rings together can produce
        // self-intersections and inconsistency with cells unaffected by the orphanage.
        // Solution: Walk the new ring and make sure its vertices are consistent with
        // the old ring.
        BGring newRing;
        std::vector<IntPoint> consistentPoints;
        for (typename BGring::const_iterator itr = thisRing.begin();
             itr != thisRing.end() - 1;
             ++itr) {

          // First check for points that are exactly equal
          bool result = false;
          typename BGring::iterator oItr = cellRings[thisIndex].begin();
          while (oItr != cellRings[thisIndex].end() and !result) {
            result += (*itr == *oItr);
            ++oItr;
          }
          if (result) consistentPoints.push_back(*(oItr-1));
          
          // No equal points. Check for points that are VERY close by. 
          else {
            RealPoint pOld, pNew = coords.dequantize(&(*itr).x);
            RealType dist = std::numeric_limits<RealType>::max();
            oItr = cellRings[thisIndex].begin();
            while (oItr != cellRings[thisIndex].end() and dist > 15.0*coords.delta) {
              pOld = coords.dequantize(&(*oItr).x);
              dist = geometry::distance<2,RealType>(&pOld.x, &pNew.x);
              ++oItr;
            }

            // No points are close by. This is likely a new ring vertex resulting
            // from the orphan, so we'll just add it.
            if (oItr == cellRings[thisIndex].end())  consistentPoints.push_back(*itr     );

            // A point is very close by. This is likely the point we were aiming for.
            // Replace the new ring vertex with the old vertex location.
            else                                     consistentPoints.push_back(*(oItr-1));
          }
        }
        consistentPoints.push_back(consistentPoints[0]);
        boost::geometry::assign(newRing, BGring(consistentPoints.begin(),
                                                consistentPoints.end()));
        boost::geometry::unique(newRing);
        thisRing = newRing;

        // // Blago!
        // std::cerr << "\nCell " << thisIndex << " after fixup" << std::endl;
        // for (typename BGring::const_iterator itr = thisRing.begin();
        //      itr != thisRing.end();
        //      ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
        // // Blago!
      }
      
      // If the orphan has only a single neighbor, just compute its union with
      // that neighbor's cell ring from the full tessellation
      else{
        thisRing = orphan;
	
	// Union orphan with its one neighbor
	std::vector<BGring> unionRing;
	boost::geometry::union_(thisRing, cellRings[thisIndex], unionRing);
	if (unionRing.size() > 1) {
	  std::cerr << "Blago!" << std::endl << "Cell " << thisIndex
		    << " has more than one cell ring:" << std::endl;
	  for (i = 0; i < unionRing.size(); ++i) {
	    std::cerr << std::endl << "Ring " << i << ":" << std::endl;
	    for (typename BGring::const_iterator itr = unionRing[i].begin();
		 itr != unionRing[i].end(); 
		 ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
	  }
	}
	POLY_ASSERT(unionRing.size() == 1);
	thisRing = unionRing[0];

#ifndef NDEBUG
        for (typename BGring::const_iterator itr = thisRing.begin();
             itr != thisRing.end() - 1; 
             ++itr ) {
          bool result = false;
          for (typename BGring::const_iterator pItr = orphan.begin();
               pItr != orphan.end() - 1;
               ++pItr) result += (*itr == *pItr);
          for (typename BGring::const_iterator pItr = cellRings[thisIndex].begin();
	       pItr != cellRings[thisIndex].end() - 1;
	       ++pItr) result += (*itr == *pItr);
          POLY_ASSERT2(result, "Union error: there are points in the single-neighbor "
                       << "union ring that are not in the orphan or its neighbor.");
        }
#endif
      }
      POLY_ASSERT(thisRing.size()  >  2);
      POLY_ASSERT(thisRing.front() == thisRing.back());
            
      // Union may produce cell rings containing edges broken into multiple
      // pieces (i.e. three or more sequential collinear nodes).
      cleanRingEdges(thisRing);
      POLY_ASSERT(thisRing.size()  >  2);
      POLY_ASSERT(thisRing.front() == thisRing.back());
      POLY_ASSERT2(!boost::geometry::intersects(thisRing), "Cell " << thisIndex);

      // // Blago!
      // std::cerr << "\nCell " << thisIndex << " after final cleaning" << std::endl;
      // for (typename BGring::const_iterator itr = thisRing.begin();
      //      itr != thisRing.end();
      //      ++itr)   std::cerr << (*itr) << coords.dequantize(&(*itr).x) << std::endl;
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
