#include "polytope_c.h"
#include "polytope.hh"

using namespace std;
using namespace polytope;

extern "C"
{

struct polytope_plc_t
{
  PLC<2, polytope_real_t>* plc2;
  PLC<3, polytope_real_t>* plc3;
};

//------------------------------------------------------------------------
polytope_plc_t* polytope_plc_new(int dimension)
{
  POLY_ASSERT(dimension == 2 or dimension == 3);
  polytope_plc_t* plc = (polytope_plc_t*)malloc(sizeof(polytope_plc_t));
  if (dimension == 2)
  {
    plc->plc2 = new PLC<2, polytope_real_t>();
    plc->plc3 = NULL;
  }
  else
  {
    plc->plc2 = NULL;
    plc->plc3 = new PLC<3, polytope_real_t>();
  }
  return plc;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_free(polytope_plc_t* plc)
{
  if (plc->plc2 != NULL)
    delete plc->plc2;
  else
    delete plc->plc3;
  free(plc);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_add_facet(polytope_plc_t* plc)
{
  vector<vector<int> >& facets = (plc->plc2 != NULL) ? plc->plc2->facets 
                                                     : plc->plc3->facets;
  facets.push_back(vector<int>());
  return facets.size()-1;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_add_facet_node(polytope_plc_t* plc, int facet, int node)
{
  size_t f = (size_t)facet;
  vector<vector<int> >& facets = (plc->plc2 != NULL) ? plc->plc2->facets 
                                                     : plc->plc3->facets;
  POLY_ASSERT(f < facets.size());
  facets[f].push_back(node);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_add_hole(polytope_plc_t* plc)
{
  vector<vector<vector<int> > >& holes = (plc->plc2 != NULL) ? plc->plc2->holes 
                                                             : plc->plc3->holes;
  holes.push_back(vector<vector<int> >());
  return holes.size()-1;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_add_hole_facet(polytope_plc_t* plc, int hole)
{
  size_t h = (size_t)hole;
  vector<vector<vector<int> > >& holes = (plc->plc2 != NULL) ? plc->plc2->holes 
                                                             : plc->plc3->holes;
  POLY_ASSERT(h < holes.size());
  holes[h].push_back(vector<int>());
  return holes[h].size()-1;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_add_hole_facet_node(polytope_plc_t* plc, int hole, int facet, int node)
{
  size_t h = (size_t)hole;
  vector<vector<vector<int> > >& holes = (plc->plc2 != NULL) ? plc->plc2->holes 
                                                             : plc->plc3->holes;
  POLY_ASSERT(h < holes.size());
  size_t f = (size_t)facet;
  POLY_ASSERT(f < holes[h].size());
  holes[h][f].push_back(node);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_num_facet_nodes(polytope_plc_t* plc, int facet)
{
  size_t f = (size_t)facet;
  vector<vector<int> >& facets = (plc->plc2 != NULL) ? plc->plc2->facets 
                                                     : plc->plc3->facets;
  POLY_ASSERT(f < facets.size());
  return (int)(facets[f].size());
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_get_facet_nodes(polytope_plc_t* plc, int facet, int* facet_nodes)
{
  size_t f = (size_t)facet;
  vector<vector<int> >& facets = (plc->plc2 != NULL) ? plc->plc2->facets 
                                                     : plc->plc3->facets;
  POLY_ASSERT(f < facets.size());
  copy(facets[f].begin(), facets[f].end(), facet_nodes);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_num_hole_facets(polytope_plc_t* plc, int hole)
{
  size_t h = (size_t)hole;
  vector<vector<vector<int> > >& holes = (plc->plc2 != NULL) ? plc->plc2->holes 
                                                             : plc->plc3->holes;
  POLY_ASSERT(h < holes.size());
  return (int)(holes[h].size());
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_num_hole_facet_nodes(polytope_plc_t* plc, int hole, int hole_facet)
{
  size_t h = (size_t)hole;
  vector<vector<vector<int> > >& holes = (plc->plc2 != NULL) ? plc->plc2->holes 
                                                             : plc->plc3->holes;
  POLY_ASSERT(h < holes.size());
  size_t f = (size_t)hole_facet;
  POLY_ASSERT(f < holes[h].size());
  return (int)(holes[h][f].size());
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_get_hole_facet_nodes(polytope_plc_t* plc, int hole, int hole_facet, int* hole_facet_nodes)
{
  size_t h = (size_t)hole;
  vector<vector<vector<int> > >& holes = (plc->plc2 != NULL) ? plc->plc2->holes 
                                                             : plc->plc3->holes;
  POLY_ASSERT(h < holes.size());
  size_t f = (size_t)hole_facet;
  POLY_ASSERT(f < holes[h].size());
  copy(holes[h][f].begin(), holes[h][f].end(), hole_facet_nodes);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_clear(polytope_plc_t* plc)
{
  if (plc->plc2 != NULL)
    plc->plc2->clear();
  else
    plc->plc3->clear();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
bool polytope_plc_empty(polytope_plc_t* plc)
{
  if (plc->plc2 != NULL)
    return plc->plc2->empty();
  else
    return plc->plc3->empty();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
bool polytope_plc_valid(polytope_plc_t* plc)
{
  if (plc->plc2 != NULL)
    return plc->plc2->valid();
  else
    return plc->plc3->valid();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_snprintf(polytope_plc_t* plc, char* str, int n)
{
  // FIXME
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_fprintf(polytope_plc_t* plc, FILE* stream)
{
  POLY_ASSERT(plc != NULL);
  // Let's assume, for the moment, that the output of our PLC scales with the 
  // number of facets.
  int n = (plc->plc2 != NULL) ? 128 * plc->plc2->facets.size()
                               : 128 * plc->plc3->facets.size();
  char* str = (char*)malloc(n * sizeof(char));
  polytope_plc_snprintf(plc, str, n);
  fprintf(stream, "%s", str);
  free(str);
}
//------------------------------------------------------------------------

}

namespace polytope
{

// Here's a helper function for constructing C PLCs from C++ ones.
template <int Dimension>
void fill_plc(const polytope::PLC<Dimension, polytope_real_t>& plc,
              polytope_plc_t* c_plc)
{
  if (Dimension == 2)
  {
    c_plc->plc2->facets = plc.facets;
    c_plc->plc2->holes = plc.holes;
  }
  else
  {
    c_plc->plc3->facets = plc.facets;
    c_plc->plc3->holes = plc.holes;
  }
}

// Here's a helper function for constructing C++ PLCs from C ones.
template <int Dimension>
void fill_plc(polytope_plc_t* c_plc,
              polytope::PLC<Dimension, polytope_real_t>& plc)
{
  if (Dimension == 2)
  {
    plc.facets = c_plc->plc2->facets;
    plc.holes = c_plc->plc2->holes;
  }
  else
  {
    plc.facets = c_plc->plc3->facets;
    plc.holes = c_plc->plc3->holes;
  }
}

// Template instantiations.
template void fill_plc(const polytope::PLC<2, polytope_real_t>&, polytope_plc_t*);
template void fill_plc(const polytope::PLC<3, polytope_real_t>&, polytope_plc_t*);

template void fill_plc(polytope_plc_t*, polytope::PLC<2, polytope_real_t>&);
template void fill_plc(polytope_plc_t*, polytope::PLC<3, polytope_real_t>&);

}

