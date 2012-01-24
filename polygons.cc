#include "polygons.hh"
#include <set>
#include <float.h>

namespace polytope
{

namespace polygons
{

using namespace std;

//-------------------------------------------------------------------
void 
traverseNodes(const Mesh<2>& mesh,
              const Cell<2>& cell, 
              vector<int>& nodes)
{
  // Start with the first face, entering both nodes.
  const Face<2>& face = *mesh.face(cell.faces[0]);
  nodes.push_back(face.nodes[0]);
  nodes.push_back(face.nodes[1]);

  // Proceed through the other faces till we're done.
  set<int> facesUsed;
  while (facesUsed.size() < cell.numFaces-1)
  {
    int numFacesUsed = facesUsed.size();
    for (int f = 1; f < cell.numFaces; ++f)
    {
      if (facesUsed.find(f) == facesUsed.end())
      {
        const Face<2>& face = *mesh.face(cell.faces[f]);

        // Does it go on the back?
        if (face.nodes[0] == nodes.back())
        {
          nodes.push_back(face.nodes[1]);
          facesUsed.insert(f);
        }
        else if (face.nodes[1] == nodes.back())
        {
          nodes.push_back(face.nodes[0]);
          facesUsed.insert(f);
        }

        // How about the front?
        else if (face.nodes[0] == nodes.front())
        {
          nodes.insert(nodes.begin(), face.nodes[1]);
          facesUsed.insert(f);
        }
        else if (face.nodes[1] == nodes.front())
        {
          nodes.insert(nodes.begin(), face.nodes[0]);
          facesUsed.insert(f);
        }
      }
    }

    // If we're here and haven't processed any more 
    // faces, we have a face remaining that doesn't 
    // connect to the end of the chain. What about 
    // the beginning?
    if (numFacesUsed == facesUsed.size())
    {
      for (int f = 1; f < cell.numFaces; ++f)
      {
        if (facesUsed.find(f) == facesUsed.end())
        {
          const Face<2>& face = *mesh.face(cell.faces[f]);
          if (face.nodes[0] == nodes.front())
          {
            nodes.push_back(face.nodes[0]);
            facesUsed.insert(f);
          }
          else if (face.nodes[1] == nodes.front())
          {
            nodes.push_back(face.nodes[1]);
            facesUsed.insert(f);
          }
        }
      }
    }
  }
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void 
traverseConvexHull(const vector<Point<2> >& points,
                   vector<int>& indices)
{
  // Find the "lowest" point in the set.
  Real ymin = FLT_MAX;
  int index0 = -1;
  for (int p = 0; p < points.size(); ++p)
  {
    if (ymin > points[p][1])
    {
      ymin = points[p][1];
      index0 = p;
    }
  }

  // We start with this point and a horizontal angle.
  Real thetaPrev = 0.0;
  indices.push_back(index0);

  // Now start gift wrapping.
  int i = index0;
  do 
  {
    Real dthetaMin = 2.0*M_PI;
    int jMin = -1;
    for (int j = 0; j < points.size(); ++j)
    {
      if (j != i)
      {
        Real dx = points[j][0] - points[i][0],
             dy = points[j][1] - points[i][1];
        Real theta = atan2(dy, dx);
        Real dtheta = theta - thetaPrev;
        if (dtheta < 0.0)
          dtheta += 2.0*M_PI;
        if (dthetaMin > dtheta)
        {
          dthetaMin = dtheta;
          jMin = j;
        }
      }
    }
    if (jMin != index0)
      indices.push_back(jMin);
    thetaPrev += dthetaMin;
    i = jMin;
  }
  while (i != index0);

  // The convex hull should be a polygon unless the input points 
  // don't form a polygon.
  ASSERT((points.size() <= 2) or 
         ((points.size() > 2) and (indices.size() > 2)));
}
//-------------------------------------------------------------------

}

}

