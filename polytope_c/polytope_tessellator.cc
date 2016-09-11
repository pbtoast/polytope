#include "polytope_c.h"
#include "polytope.hh"
#include "BoostTessellator.hh"
#include "TriangleTessellator.hh"
#include "TetgenTessellator.hh"

namespace polytope
{

// Helper functions for constructing C PLCs from C++ ones and vice versa. 
// Defined in polytope_plc.cc.
template <int Dimension>
void fill_plc(const polytope::PLC<Dimension, polytope_real_t>& plc,
              polytope_plc_t* c_plc);

template <int Dimension>
void fill_plc(polytope_plc_t* c_plc,
              polytope::PLC<Dimension, polytope_real_t>& plc);

// This helper crafts a C polytope_tessellation object from a C++ Tessellation.
template <int Dimension>
void fill_tessellation(const Tessellation<Dimension, polytope_real_t>& t, polytope_tessellation_t* tess)
{
  // Copy node coordinates.
  tess->num_nodes = t.nodes.size()/Dimension;
  tess->nodes = (polytope_real_t*)malloc(sizeof(polytope_real_t) * Dimension * tess->num_nodes);
  copy(t.nodes.begin(), t.nodes.end(), tess->nodes);

  // Copy cell-face data.
  {
    tess->num_cells = t.cells.size();
    tess->cell_offsets = (int*)malloc(sizeof(int) * (tess->num_cells+1));
    int size = 0;
    for (size_t i = 0; i < t.cells.size(); ++i)
      size += t.cells[i].size();
    tess->cell_faces = (int*)malloc(sizeof(int) * (size+1));
    int offset = 0;
    for (size_t i = 0; i < t.cells.size(); ++i)
    {
      tess->cell_offsets[i] = offset;
      for (size_t j = 0; j < t.cells[i].size(); ++j, ++offset)
        tess->cell_faces[offset] = t.cells[i][j];
    }
    tess->cell_offsets[t.cells.size()] = offset;
  }

  // Copy face-node data.
  {
    tess->num_faces = t.faces.size();
    tess->face_offsets = (int*)malloc(sizeof(int) * (tess->num_faces+1));
    int size = 0;
    for (size_t i = 0; i < t.faces.size(); ++i)
      size += t.faces[i].size();
    tess->face_nodes = (unsigned*)malloc(sizeof(unsigned) * (size+1));
    int offset = 0;
    for (size_t i = 0; i < t.faces.size(); ++i)
    {
      tess->face_offsets[i] = offset;
      for (size_t j = 0; j < t.faces[i].size(); ++j, ++offset)
        tess->face_nodes[offset] = t.faces[i][j];
    }
    tess->face_offsets[t.faces.size()] = offset;
  }

  // Boundary nodes and faces.
  tess->num_boundary_nodes = (int)t.boundaryNodes.size();
  tess->boundary_nodes = (unsigned*)malloc(sizeof(unsigned) * tess->num_boundary_nodes);
  copy(t.boundaryNodes.begin(), t.boundaryNodes.end(), tess->boundary_nodes);

  tess->num_boundary_faces = (int)t.boundaryFaces.size();
  tess->boundary_faces = (unsigned*)malloc(sizeof(unsigned) * tess->num_boundary_faces);
  copy(t.boundaryFaces.begin(), t.boundaryFaces.end(), tess->boundary_faces);

  // Cells attached to faces.
  POLY_ASSERT(tess->num_faces == (int)t.faceCells.size());
  tess->face_cells = (int*)malloc(2 * tess->num_faces * sizeof(int));
  for (int f = 0; f < t.faceCells.size(); ++f)
  {
    tess->face_cells[2*f] = t.faceCells[f][0];
    tess->face_cells[2*f+1] = (t.faceCells[f].size() == 2) ? t.faceCells[f][1] : -1;
  }

  // Convex hull.
  if (!t.convexHull.empty())
  {
    tess->convex_hull = polytope_plc_new(Dimension);
    fill_plc(t.convexHull, tess->convex_hull);
  }

  // Neighbor domain information.
  tess->num_neighbor_domains = (int)t.neighborDomains.size();
  tess->neighbor_domains = (unsigned*)malloc(sizeof(unsigned) * tess->num_neighbor_domains);
  copy(t.neighborDomains.begin(), t.neighborDomains.end(), tess->neighbor_domains);

  {
    tess->shared_node_domain_offsets = (int*)malloc(sizeof(int) * (tess->num_neighbor_domains + 1));
    int size = 0;
    for (size_t i = 0; i < t.sharedNodes.size(); ++i)
      size += t.sharedNodes[i].size();
    tess->shared_nodes = (unsigned*)malloc(sizeof(unsigned) * (size+1));
    int offset = 0;
    for (size_t i = 0; i < t.sharedNodes.size(); ++i)
    {
      tess->shared_node_domain_offsets[i] = offset;
      for (size_t j = 0; j < t.sharedNodes[i].size(); ++j, ++offset)
        tess->shared_nodes[offset] = t.sharedNodes[i][j];
    }
    tess->shared_node_domain_offsets[t.sharedNodes.size()] = offset;
  }

  {
    tess->shared_face_domain_offsets = (int*)malloc(sizeof(int) * (tess->num_neighbor_domains + 1));
    int size = 0;
    for (size_t i = 0; i < t.sharedFaces.size(); ++i)
      size += t.sharedFaces[i].size();
    tess->shared_faces = (unsigned*)malloc(sizeof(unsigned) * (size+1));
    int offset = 0;
    for (size_t i = 0; i < t.sharedFaces.size(); ++i)
    {
      tess->shared_face_domain_offsets[i] = offset;
      for (size_t j = 0; j < t.sharedFaces[i].size(); ++j, ++offset)
        tess->shared_faces[offset] = t.sharedFaces[i][j];
    }
    tess->shared_face_domain_offsets[t.sharedFaces.size()] = offset;
  }


  // Node->cell connectivity.
  // FIXME
}

// This helper crafts a C++ Tessellation object from a C polytope_tessellation_t.
template <int Dimension>
void fill_tessellation(polytope_tessellation_t* tess, Tessellation<Dimension, polytope_real_t>& t)
{
  // Copy node coordinates.
  t.nodes.resize(Dimension * tess->num_nodes);
  copy(tess->nodes, tess->nodes + Dimension * tess->num_nodes, t.nodes.begin());

  // Copy cell-face data.
  t.cells.resize(tess->num_cells);
  for (size_t i = 0; i < t.cells.size(); ++i)
  {
    int num_faces = tess->cell_offsets[i+1] - tess->cell_offsets[i];
    t.cells[i].resize(num_faces);
    for (int j = 0; j < num_faces; ++j)
      t.cells[i][j] = tess->cell_faces[tess->cell_offsets[i] + j];
  }

  // Copy face-node data.
  t.faces.resize(tess->num_faces);
  for (size_t i = 0; i < t.faces.size(); ++i)
  {
    int num_nodes = tess->face_offsets[i+1] - tess->face_offsets[i];
    t.faces[i].resize(num_nodes);
    for (int j = 0; j < num_nodes; ++j)
      t.faces[i][j] = tess->face_nodes[tess->face_offsets[i] + j];
  }

  // Boundary nodes and faces.
  t.boundaryNodes.resize(tess->num_boundary_nodes);
  copy(tess->boundary_nodes, tess->boundary_nodes + tess->num_boundary_nodes, t.boundaryNodes.begin());

  t.boundaryFaces.resize(tess->num_boundary_faces);
  copy(tess->boundary_faces, tess->boundary_faces + tess->num_boundary_faces, t.boundaryFaces.begin());

  // Cells attached to faces.
  t.faceCells.resize(tess->num_faces);
  for (size_t i = 0; i < t.faceCells.size(); ++i)
  {
    int num_cells = (tess->face_cells[2*i+1] != -1) ? 2 : 1;
    t.faceCells[i].resize(num_cells);
    for (int j = 0; j < num_cells; ++j)
      t.faceCells[i][j] = tess->face_cells[2*i+j];
  }

  // Convex hull.
  // FIXME: Not yet supported.

  // Neighbor domain information.
  t.neighborDomains.resize(tess->num_neighbor_domains);
  t.sharedNodes.resize(tess->num_neighbor_domains);
  t.sharedFaces.resize(tess->num_neighbor_domains);
  copy(tess->neighbor_domains, tess->neighbor_domains + tess->num_neighbor_domains, t.neighborDomains.begin());
  for (size_t i = 0; i < t.neighborDomains.size(); ++i)
  {
    int num_nodes = tess->shared_node_domain_offsets[i+1] - tess->shared_node_domain_offsets[i];
    t.sharedNodes[i].resize(num_nodes);
    for (int j = 0; j < num_nodes; ++j)
      t.sharedNodes[i][j] = tess->shared_nodes[tess->shared_node_domain_offsets[i] + j];

    int num_faces = tess->shared_face_domain_offsets[i+1] - tess->shared_face_domain_offsets[i];
    t.sharedFaces[i].resize(num_faces);
    for (int j = 0; j < num_faces; ++j)
      t.sharedFaces[i][j] = tess->shared_faces[tess->shared_face_domain_offsets[i] + j];
  }

  // We ignore Node->cell connectivity, since it can always be generated.
}

// Template instantiations.
template void fill_tessellation(const Tessellation<2, polytope_real_t>& t, polytope_tessellation_t* tess);
template void fill_tessellation(const Tessellation<3, polytope_real_t>& t, polytope_tessellation_t* tess);
template void fill_tessellation(polytope_tessellation_t* tess, Tessellation<2, polytope_real_t>& t);
template void fill_tessellation(polytope_tessellation_t* tess, Tessellation<3, polytope_real_t>& t);

}

using namespace std;
using namespace polytope;

extern "C"
{

//------------------------------------------------------------------------
struct polytope_tessellator_t
{
  Tessellator<2, polytope_real_t>* tess2;
  Tessellator<3, polytope_real_t>* tess3;
};
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellator_free(polytope_tessellator_t* tessellator)
{
  if (tessellator->tess2 != NULL)
    delete tessellator->tess2;
  else
    delete tessellator->tess3;
  free(tessellator);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellator_tessellate_unbounded(polytope_tessellator_t* tessellator,
                                               polytope_real_t* points, int num_points,
                                               polytope_tessellation_t* mesh)
{
  if (tessellator->tess2 != NULL)
  {
    vector<polytope_real_t> pts(points, points + 2*num_points);
    Tessellation<2, polytope_real_t> t;
    tessellator->tess2->tessellate(pts, t);
    fill_tessellation(t, mesh);
  }
  else
  {
    POLY_ASSERT(tessellator->tess3 != NULL);
    vector<polytope_real_t> pts(points, points + 3*num_points);
    Tessellation<3, polytope_real_t> t;
    tessellator->tess3->tessellate(pts, t);
    fill_tessellation(t, mesh);
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellator_tessellate_in_box(polytope_tessellator_t* tessellator,
                                            polytope_real_t* points, int num_points,
                                            polytope_real_t* low,
                                            polytope_real_t* high,
                                            polytope_tessellation_t* mesh)
{
  if (tessellator->tess2 != NULL)
  {
    vector<polytope_real_t> pts(points, points + 2*num_points);
    Tessellation<2, polytope_real_t> t;
    tessellator->tess2->tessellate(pts, low, high, t);
    fill_tessellation(t, mesh);
  }
  else
  {
    POLY_ASSERT(tessellator->tess3 != NULL);
    vector<polytope_real_t> pts(points, points + 3*num_points);
    Tessellation<3, polytope_real_t> t;
    tessellator->tess3->tessellate(pts, low, high, t);
    fill_tessellation(t, mesh);
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellator_tessellate_in_plc(polytope_tessellator_t* tessellator,
                                            polytope_real_t* points, int num_points,
                                            polytope_real_t* plc_points, int num_plc_points,
                                            polytope_plc_t* piecewise_linear_complex,
                                            polytope_tessellation_t* mesh)
{
  if (tessellator->tess2 != NULL)
  {
    vector<polytope_real_t> pts(points, points + 2*num_points);
    vector<polytope_real_t> plc_pts(plc_points, plc_points + 2*num_plc_points);
    PLC<2, polytope_real_t> plc;
    fill_plc(piecewise_linear_complex, plc);
    Tessellation<2, polytope_real_t> t;
    tessellator->tess2->tessellate(pts, plc_pts, plc, t);
    fill_tessellation(t, mesh);
  }
  else
  {
    POLY_ASSERT(tessellator->tess3 != NULL);
    vector<polytope_real_t> pts(points, points + 3*num_points);
    vector<polytope_real_t> plc_pts(plc_points, plc_points + 3*num_plc_points);
    PLC<3, polytope_real_t> plc;
    fill_plc(piecewise_linear_complex, plc);
    Tessellation<3, polytope_real_t> t;
    tessellator->tess3->tessellate(pts, plc_pts, plc, t);
    fill_tessellation(t, mesh);
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
bool polytope_tessellator_handles_plcs(polytope_tessellator_t* tessellator)
{
  if (tessellator->tess2 != NULL)
    return tessellator->tess2->handlesPLCs();
  else
    return tessellator->tess3->handlesPLCs();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
const char* polytope_tessellator_name(polytope_tessellator_t* tessellator)
{
  if (tessellator->tess2 != NULL)
    return tessellator->tess2->name().c_str();
  else
    return tessellator->tess3->name().c_str();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
polytope_real_t polytope_tessellator_degeneracy(polytope_tessellator_t* tessellator)
{
  if (tessellator->tess2 != NULL)
    return tessellator->tess2->degeneracy();
  else
    return tessellator->tess3->degeneracy();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_tessellator_dimension(polytope_tessellator_t* tessellator)
{
  if (tessellator->tess2 != NULL)
    return 2;
  else
    return 3;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
#ifdef HAVE_BOOST
polytope_tessellator_t* boost_tessellator_new()
{
  polytope_tessellator_t* t = (polytope_tessellator_t*)malloc(sizeof(polytope_tessellator_t));
  t->tess2 = new BoostTessellator<polytope_real_t>();
  t->tess3 = NULL;
  return t;
}
#endif
//------------------------------------------------------------------------

//------------------------------------------------------------------------
#ifdef HAVE_TRIANGLE
polytope_tessellator_t* triangle_tessellator_new()
{
  polytope_tessellator_t* t = (polytope_tessellator_t*)malloc(sizeof(polytope_tessellator_t));
  t->tess2 = new TriangleTessellator<polytope_real_t>();
  t->tess3 = NULL;
  return t;
}
#endif
//------------------------------------------------------------------------

//------------------------------------------------------------------------
#ifdef HAVE_TETGEN
polytope_tessellator_t* tetgen_tessellator_new()
{
  polytope_tessellator_t* t = (polytope_tessellator_t*)malloc(sizeof(polytope_tessellator_t));
  t->tess2 = NULL;
  t->tess3 = new TetgenTessellator();
  return t;
}
#endif
//------------------------------------------------------------------------

}

