//------------------------------------------------------------------------
// TriangleTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <set>
#include "float.h"
#include <iostream>

#include "polytope.hh" // Pulls in ASSERT and TriangleTessellator.hh.

// Since triangle isn't built to work out-of-the-box with C++, we 
// slurp in its source here, bracketing it with the necessary dressing.
#define TRILIBRARY
#define REAL double
#define ANSI_DECLARATORS 
#define CDT_ONLY // Conforming Delaunay triangulations only! 

// Because Hang Si has messed with some of Shewchuk's predicates and
// included them with his own Tetgen library, we need to rename some of 
// the symbols therein to prevent duplicate symbols from confusing the 
// linker. Gross.
#define exactinit triangle_exactinit
#define fast_expansion_sum_zeroelim triangle_fast_expansion_sum_zeroelim
#define scale_expansion_zeroelim triangle_scale_expansion_zeroelim
#define estimate triangle_estimate
#define orient3dadapt triangle_orient3dadapt
#define orient3d triangle_orient3d
#define incircleadapt triangle_incircleadapt
#define incircle triangle_incircle
#include "triangle.c"

// Fast predicate for determining colinearity of points.
extern double orient2d(double* pa, double* pb, double* pc);

namespace polytope {

using namespace std;

namespace {

//------------------------------------------------------------------------
// This function computes the circumcenter of a triangle with vertices
// A = (Ax, Ay), B = (Bx, By), and C = (Cx, Cy), and places the result 
// in X.
//------------------------------------------------------------------------
void 
computeCircumcenter(double* A, double* B, double* C, double* X)
{
  // This solution was taken from Wikipedia's entry:
  // http://en.wikipedia.org/wiki/Circumscribed_circle
  double D = 2.0*(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1]));
  X[0] = ((A[0]*A[0] + A[1]*A[1])*(B[1]-C[1]) + (B[0]*B[0] + B[1]*B[1])*(C[1]-A[1]) + 
          (C[0]*C[0] + C[1]*C[1])*(A[1]-B[1]))/D;
  X[1] = ((A[0]*A[0] + A[1]*A[1])*(C[0]-B[0]) + (B[0]*B[0] + B[1]*B[1])*(A[0]-C[0]) + 
          (C[0]*C[0] + C[1]*C[1])*(B[0]-A[0]))/D;
}
//------------------------------------------------------------------------

} // end anonymous namespace

//------------------------------------------------------------------------------
template<typename Real>
TriangleTessellator<Real>::
TriangleTessellator():
  Tessellator<2, Real>()
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
TriangleTessellator<Real>::
~TriangleTessellator() 
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
void
TriangleTessellator<Real>::
tessellate(const vector<Real>& points,
           Tessellation<2, Real>& mesh) const 
{
  // Create an empty PLC and go for it.
  PLC<2, Real> noBoundary;
  tessellate(points, noBoundary, mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Real>
void
TriangleTessellator<Real>::
tessellate(const vector<Real>& points,
           const PLC<2, Real>& geometry,
           Tessellation<2, Real>& mesh) const 
{
  ASSERT(!points.empty());

  // Make sure we're not modifying an existing tessellation.
  ASSERT(mesh.empty());

  triangulateio in, delaunay;

  // Define input points.
  in.numberofpoints = points.size()/2;
  in.pointlist = new Real[points.size()];
  copy(points.begin(), points.end(), in.pointlist);

  // No point attributes or markers.
  in.numberofpointattributes = 0;
  in.pointattributelist = 0; 
  in.pointmarkerlist = 0; // No point markers.

  // Segments and/or holes.
  if (geometry.empty())
  {
    in.numberofsegments = 0;
    in.segmentlist = 0;
    in.segmentmarkerlist = 0;
    in.numberofholes = 0;
    in.holelist = 0;
  }
  else
  {
    in.numberofsegments = geometry.facets.size();
    in.segmentlist = new int[2*in.numberofsegments];
    int s = 0;
    for (int f = 0; f < geometry.facets.size(); ++f)
    {
      in.segmentlist[s++] = geometry.facets[f][0];
      in.segmentlist[s++] = geometry.facets[f][1];
    }
    in.segmentmarkerlist = 0;
    in.numberofholes = geometry.holes.size();
    in.holelist = new double[2*in.numberofholes];
    copy(geometry.holes.begin(), geometry.holes.end(), in.holelist);
  }

  // No regions.
  in.numberofregions = 0;
  in.regionlist = 0;

  // Set up the structure for the triangulation.
  delaunay.pointlist = 0;
  delaunay.pointattributelist = 0;
  delaunay.pointmarkerlist = 0;
  delaunay.trianglelist = 0;
  delaunay.triangleattributelist = 0;
  delaunay.neighborlist = 0;
  if (geometry.empty())
  {
    delaunay.segmentlist = 0;
    delaunay.segmentmarkerlist = 0;
  }
  else
  {
    delaunay.numberofsegments = geometry.facets.size();
    delaunay.segmentlist = new int[2*delaunay.numberofsegments];
    delaunay.segmentmarkerlist = new int[delaunay.numberofsegments];
  }
  delaunay.edgelist = 0;
  delaunay.edgemarkerlist = 0;
  delaunay.holelist = 0;

  // Do the triangulation. Switches pass to triangle are:
  // -Q : Quiet (shaddap!), no output on the terminal except errors.
  // -z : Indices are all numbered from zero.
  // -e : Generates edges and places them in out.edgelist.
  // -c : Generates convex hull and places it in out.segmentlist.
  // -p : Uses the given PLC information.
  if (geometry.empty())
    triangulate((char*)"Qzec", &in, &delaunay, 0);
  else
    triangulate((char*)"Qzep", &in, &delaunay, 0);

  // Make sure we got something.
  ASSERT(delaunay.numberoftriangles > 0);

  //--------------------------------------------------------
  // Create the Voronoi tessellation from the triangulation.
  //--------------------------------------------------------

  // Create the Voronoi nodes from the list of triangles. Each triangle 
  // has 3 nodes p, q, r, and corresponds to a Voronoi node at (X,Y), say.
  // The Voronoi node is located at the center of the circle that contains 
  // p, q, and r. 
  vector<vector<int> > cellNodes(delaunay.numberofpoints);

  // A table of Voronoi nodes falling on Delaunay edges. This prevents
  // us from double-counting nodes.
  map<pair<int, int>, int> nodesOnEdges; 
  for (int i = 0; i < delaunay.numberoftriangles; ++i)
  {
    // Coordinates of the triangle's vertices p, q, r.
    int pindex = delaunay.trianglelist[3*i];
    double p[2];
    p[0] = delaunay.pointlist[2*pindex];
    p[1] = delaunay.pointlist[2*pindex+1];

    int qindex = delaunay.trianglelist[3*i+1];
    double q[2];
    q[0] = delaunay.pointlist[2*qindex],
    q[1] = delaunay.pointlist[2*qindex+1];

    int rindex = delaunay.trianglelist[3*i+2];
    double r[2];
    r[0] = delaunay.pointlist[2*rindex],
    r[1] = delaunay.pointlist[2*rindex+1];

    // Get the circumcenter.
    double X[2];
    computeCircumcenter(p, q, r, X);

    // If this node lies on an edge, add it to the table.
    // This uses Jonathan Shewchuck's fast geometry predicate orient2d().
    if (orient2d(p, q, X) == 0.0)
    {
      pair<int, int> key(min(pindex, qindex), max(pindex, qindex));
      map<pair<int, int>, int>::const_iterator iter = nodesOnEdges.find(key);
      if (iter != nodesOnEdges.end())
      {
        // Add this existing node to the missing cell and proceed.
        cellNodes[rindex].push_back(iter->second);
        continue;
      }
      else
        nodesOnEdges[key] = i;
    }
    else if (orient2d(p, r, X) == 0.0) 
    {
      pair<int, int> key(min(pindex, rindex), max(pindex, rindex));
      map<pair<int, int>, int>::const_iterator iter = nodesOnEdges.find(key);
      if (iter != nodesOnEdges.end())
      {
        // Add this existing node to the missing cell and proceed.
        cellNodes[qindex].push_back(iter->second);
        continue;
      }
      else
        nodesOnEdges[key] = i;
    }
    else if (orient2d(q, r, X) == 0.0)
    {
      pair<int, int> key(min(qindex, rindex), max(qindex, rindex));
      map<pair<int, int>, int>::const_iterator iter = nodesOnEdges.find(key);
      if (iter != nodesOnEdges.end())
      {
        // Add this existing node to the missing cell and proceed.
        cellNodes[pindex].push_back(iter->second);
        continue;
      }
      else
        nodesOnEdges[key] = i;
    }

    // Assign the node coordinate.
    mesh.nodes.push_back(Real(X[0]));
    mesh.nodes.push_back(Real(X[1]));

    // Associate this node with its Voronoi cells.
    cellNodes[pindex].push_back(i);
    cellNodes[qindex].push_back(i);
    cellNodes[rindex].push_back(i);
  }

  // Compute the Voronoi faces between cells. These have a one-to-one
  // correspondence with Delaunay edges of triangles unless we have placed 
  // nodes along edges.
  mesh.cells.resize(delaunay.numberofpoints);
  for (int i = 0; i < delaunay.numberofedges; ++i)
  {
    int cell1 = delaunay.edgelist[2*i],
        cell2 = delaunay.edgelist[2*i+1];

    // Have we placed any nodes along this edge? If so, there's no 
    // face here.
    pair<int, int> key(min(cell1, cell2), max(cell1, cell2));
    map<pair<int, int>, int>::const_iterator iter = nodesOnEdges.find(key);
    if (iter != nodesOnEdges.end())
//{
//cout << "cells " << cell1 << ", " << cell2 << " have a node on their edge." << endl;
      continue;
//}

    // Hook the cells up to the faces.
    mesh.faces.push_back(vector<unsigned>());
    mesh.faceCells.push_back(vector<unsigned>());
    mesh.faceCells.back().reserve(2);
    mesh.faceCells.back().push_back(cell1);
    mesh.faceCells.back().push_back(cell2);

    // Hook the faces up to the cells.
    mesh.cells[cell1].push_back(mesh.faces.size()-1);
    mesh.cells[cell2].push_back(mesh.faces.size()-1);

    // Nodes that are attached to both cell 1 and cell2 are nodes 
    // of this face.
    set<int> cell1Nodes(cellNodes[cell1].begin(), cellNodes[cell1].end()),
             cell2Nodes(cellNodes[cell2].begin(), cellNodes[cell2].end()),
             faceNodes;
    set_intersection(cell1Nodes.begin(), cell1Nodes.end(),
                     cell2Nodes.begin(), cell2Nodes.end(),
                     inserter(faceNodes, faceNodes.begin()));
    mesh.faces.back().resize(faceNodes.size());
//cout << "Face " << mesh.faces.size()-1 << " (between " << cell1 << " and " << cell2 << ") has " << faceNodes.size() << " nodes\n";
    copy(faceNodes.begin(), faceNodes.end(), mesh.faces.back().begin());
  }

  // If no boundary was specified, compute the convex hull and leave.
  if (geometry.empty()) 
  {
    set<unsigned> convexHull;
    for (int i = 0; i < delaunay.numberofsegments; ++i)
    {
      int cell1 = delaunay.segmentlist[2*i];
      int cell2 = delaunay.segmentlist[2*i+1];
      convexHull.insert(cell1);
      convexHull.insert(cell2);
    }
    mesh.convexHull.resize(convexHull.size());
    copy(convexHull.begin(), convexHull.end(), mesh.convexHull.begin());

    // Clean up.
    delete [] in.pointlist;
    trifree((VOID*)delaunay.pointlist);
    trifree((VOID*)delaunay.trianglelist);
    trifree((VOID*)delaunay.edgelist);
    if (geometry.empty())
    {
      trifree((VOID*)delaunay.segmentlist);
      trifree((VOID*)delaunay.segmentmarkerlist);
    }
    else
    {
      delete [] in.segmentlist;
      delete [] in.holelist;
      delete [] delaunay.segmentlist;
      delete [] delaunay.segmentmarkerlist;
    }
    return;
  }

  // At this point, all of the interior cells are set up properly, and all 
  // nodes exist. However, we still have to finish constructing the 
  // Voronoi cells that sit at the boundary, since they don't have "back" 
  // faces. This Tessellator constructs faces that coincide with the 
  // convex hull of the triangulation. The convex hull is expressed in 
  // terms of "segments" connecting generator points. 
  
  // Extract the convex hull here, and make sure each existing face 
  // has two nodes.
  set<unsigned> convexHull;
  set<pair<int, int> > patchedFaces;
  for (int i = 0; i < delaunay.numberofsegments; ++i)
  {
    int cell1 = delaunay.segmentlist[2*i];
    int cell2 = delaunay.segmentlist[2*i+1];
    convexHull.insert(cell1);
    convexHull.insert(cell2);
//cout << "Inspecting segment for " << cell1 << ", " << cell2 << endl;

    // Add nodes to all faces attached to these cell with only 
    // one node.
    for (size_t f = 0; f < mesh.cells[cell1].size(); ++f)
    {
      int face = mesh.cells[cell1][f];
//cout << "face " << face << " has " << mesh.faces[face].size() << " nodes\n";
      if (mesh.faces[face].size() == 1)
      {
        // We've found a face with a single node. Find out where the
        // other node should go and create it.
        for (size_t f1 = 0; f1 < mesh.cells[cell2].size(); ++f1)
        {
          int face1 = mesh.cells[cell2][f1];
          if (face1 == face)
          {
            // The node should bisect the segment on the convex hull 
            // between cell1 and cell2.
            Real nodex = 0.5*(points[2*cell1]+points[2*cell2]),
                 nodey = 0.5*(points[2*cell1+1]+points[2*cell2+1]);
//cout << "Adding node (" << nodex << ", " << nodey << ") for cells " << cell1 << ", " << cell2 << endl;
            mesh.nodes.push_back(nodex);
            mesh.nodes.push_back(nodey);
            mesh.faces[face].push_back(mesh.nodes.size()/2 - 1);
            break;
          }
        }
      }
    }
  }
  mesh.convexHull.resize(convexHull.size());
  copy(convexHull.begin(), convexHull.end(), mesh.convexHull.begin());

  // Create extra nodes for the boundary faces. These nodes coincide with 
  // the Voronoi generater points.
  size_t oldNumNodes = mesh.nodes.size() / 2;
//cout << "Mesh had " << oldNumNodes << ", now adding " << mesh.convexHull.size() << endl;
  mesh.nodes.resize(2*(oldNumNodes + mesh.convexHull.size()));
  for (size_t i = 0; i < mesh.convexHull.size(); ++i)
  {
    unsigned pindex = mesh.convexHull[i];
    mesh.nodes[2*(oldNumNodes+i)]   = points[2*pindex];
    mesh.nodes[2*(oldNumNodes+i)+1] = points[2*pindex+1];
//cout << "Adding node (" << points[2*pindex] << ", " << points[2*pindex+1] << ")\n";
  }

  // Now we construct the remaining faces for the boundary cells.
  for (int i = 0; i < mesh.convexHull.size(); ++i)
  {
    int cell = mesh.convexHull[i]; // The boundary cell.

    // Add two new faces for this cell.
    mesh.faces.resize(mesh.faces.size()+2);
    mesh.faceCells.resize(mesh.faceCells.size()+2);

    // Each face has one node that it shares with an existing face.
    // There are two such existing faces, and each has a node that 
    // it does not yet share with any other.
    map<unsigned, int> numFacesForNode;
    for (size_t f = 0; f < mesh.cells[cell].size(); ++f)
    {
      int face = mesh.cells[cell][f];
      for (size_t n = 0; n < mesh.faces[face].size(); ++n)
        ++numFacesForNode[mesh.faces[face][n]];
    }
    int node1 = -1, node2 = -1;
    for (map<unsigned, int>::const_iterator iter = numFacesForNode.begin();
         iter != numFacesForNode.end(); ++iter)
    {
      if ((iter->second < 2) and ((node1 == -1) or (node2 == -1)))
      {
        if (node1 == -1)
        {
          node1 = iter->first;
        }
        else
        {
          ASSERT(node2 == -1);
          node2 = iter->first;
        }
      }
    }
    ASSERT(node1 >= 0);
    ASSERT(node2 >= 0);
    mesh.faces[mesh.faces.size()-2].push_back(node1);
    mesh.faces[mesh.faces.size()-1].push_back(node2);

    // The remaining node on each face is the one coinciding with 
    // the generator point.
    int genPointNode = oldNumNodes+i;
    mesh.faces[mesh.faces.size()-2].push_back(genPointNode);
    mesh.faces[mesh.faces.size()-1].push_back(genPointNode);

    // Each face has only the 1 cell.
    mesh.faceCells[mesh.faceCells.size()-2].push_back(cell);
    mesh.faceCells[mesh.faceCells.size()-1].push_back(cell);
    mesh.cells[cell].push_back(mesh.faces.size()-2);
    mesh.cells[cell].push_back(mesh.faces.size()-1);
  }

  // Clean up.
  delete [] in.pointlist;
  trifree((VOID*)delaunay.pointlist);
  trifree((VOID*)delaunay.trianglelist);
  trifree((VOID*)delaunay.edgelist);
  if (geometry.empty())
  {
    trifree((VOID*)delaunay.segmentlist);
    trifree((VOID*)delaunay.segmentmarkerlist);
  }
  else
  {
    delete [] in.segmentlist;
    delete [] in.holelist;
    delete [] delaunay.segmentlist;
    delete [] delaunay.segmentmarkerlist;
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TriangleTessellator<double>;

}
