#include "polytope_plc.h"
#include "polytope.hh"

using namespace std;
using namespace polytope;

extern "C"
{

struct polytope_plc_t
{
  int dimension;
  PLC<2, polytope_real_t> plc2;
  PLC<3, polytope_real_t> plc3;
}

//------------------------------------------------------------------------
polytope_plc_t* polytope_plc_new(int dimension)
{
  POLY_ASSERT(dimension == 2 or dimension == 3);
  polytope_plc_t* plc = (polytope_plc_t*)malloc(sizeof(polytope_plc_t));
  plc->dimension = dimension;
  return plc;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_free(polytope_plc_t* plc)
{
  free(plc);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_add_facet(polytope_plc_t* plc)
{
  if (plc->dimension == 2)
  {
    plc->plc2.facets.push_back(vector<int>());
    return plc->plc2.facets.size()-1;
  }
  else
  {
    plc->plc3.facets.push_back(vector<int>());
    return plc->plc3.facets.size()-1;
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_add_facet_node(polytope_plc_t* plc, int facet, int node)
{
  if (plc->dimension == 2)
  {
    size_t f = (size_t)facet;
    POLY_ASSERT(f < plc->plc2.facets.size());
    plc->plc2.facets[f].push_back(node);
  }
  else
  {
    size_t f = (size_t)facet;
    POLY_ASSERT(f < plc->plc3.facets.size());
    plc->plc3.facets[f].push_back(node);
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_add_hole(polytope_plc_t* plc)
{
  if (plc->dimension == 2)
  {
    plc->plc2.holes.push_back(vector<int>());
    return plc->plc2.holes.size()-1;
  }
  else
  {
    plc->plc3.holes.push_back(vector<int>());
    return plc->plc3.holes.size()-1;
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_add_hole_facet(polytope_plc_t* plc, int hole)
{
  if (plc->dimension == 2)
  {
    size_t h = (size_t)hole;
    POLY_ASSERT(h < plc->plc2.holes.size());
    plc->plc2.holes[h].push_back(vector<int>());
    return plc->plc2.holes[h].size()-1;
  }
  else
  {
    size_t f = (size_t)facet;
    POLY_ASSERT(h < plc->plc3.holes.size());
    plc->plc3.holes[h].push_back(vector<int>());
    return plc->plc3.holes[h].size()-1;
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_add_hole_facet_node(polytope_plc_t* plc, int facet, int node)
{
  if (plc->dimension == 2)
  {
    size_t h = (size_t)hole;
    POLY_ASSERT(h < plc->plc2.holes.size());
    size_t f = (size_t)facet;
    POLY_ASSERT(f < plc->plc2.holes[h].size());
    plc->plc2.holes[h][f].push_back(node);
  }
  else
  {
    size_t h = (size_t)hole;
    POLY_ASSERT(h < plc->plc3.holes.size());
    size_t f = (size_t)facet;
    POLY_ASSERT(f < plc->plc3.holes[h].size());
    plc->plc3.holes[h][f].push_back(node);
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_num_facet_nodes(polytope_plc_t* plc, int facet)
{
  if (plc->dimension == 2)
  {
    size_t f = (size_t)facet;
    POLY_ASSERT(f < plc->plc2.facets.size());
    return (int)plc->plc2.facets[f].size();
  }
  else
  {
    size_t f = (size_t)facet;
    POLY_ASSERT(f < plc->plc3.facets.size());
    return (int)plc->plc3.facets[f].size();
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_get_facet_nodes(polytope_plc_t* plc, int facet, int* facet_nodes)
{
  if (plc->dimension == 2)
  {
    size_t f = (size_t)facet;
    POLY_ASSERT(f < plc->plc2.facets.size());
    copy(plc->plc2.facets.begin(), plc->plc2.facets.end(), facet_nodes);
  }
  else
  {
    size_t f = (size_t)facet;
    POLY_ASSERT(f < plc->plc3.facets.size());
    copy(plc->plc3.facets.begin(), plc->plc3.facets.end(), facet_nodes);
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_num_hole_facets(polytope_plc_t* plc, int hole)
{
  if (plc->dimension == 2)
  {
    size_t h = (size_t)hole;
    POLY_ASSERT(h < plc->plc2.holes.size());
    return (int)plc->plc2.holes[h].size();
  }
  else
  {
    size_t h = (size_t)hole;
    POLY_ASSERT(h < plc->plc3.holes.size());
    return (int)plc->plc3.holes[h].size();
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_num_hole_facet_nodes(polytope_plc_t* plc, int hole, int hole_facet)
{
  if (plc->dimension == 2)
  {
    size_t h = (size_t)hole;
    POLY_ASSERT(h < plc->plc2.holes.size());
    size_t f = (size_t)hole_facet;
    POLY_ASSERT(f < plc->plc2.holes[h].size());
    return (int)plc->plc2.holes[h][f].size();
  }
  else
  {
    size_t h = (size_t)hole;
    POLY_ASSERT(h < plc->plc3.holes.size());
    size_t f = (size_t)hole_facet;
    POLY_ASSERT(f < plc->plc3.holes[h].size());
    return (int)plc->plc3.holes[h][f].size();
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_get_hole_facet_nodes(polytope_plc_t* plc, int hole, int hole_facet, int* hole_facet_nodes)
{
  if (plc->dimension == 2)
  {
    size_t h = (size_t)hole;
    POLY_ASSERT(h < plc->plc2.holes.size());
    size_t f = (size_t)hole_facet;
    POLY_ASSERT(f < plc->plc2.holes[h].size());
    copy(plc->plc2.holes[h][f].begin(), plc->plc2.holes[h][f].end(), hole_facet_nodes);
  }
  else
  {
    size_t h = (size_t)hole;
    POLY_ASSERT(h < plc->plc3.holes.size());
    size_t f = (size_t)hole_facet;
    POLY_ASSERT(f < plc->plc3.holes[h].size());
    copy(plc->plc3.holes[h][f].begin(), plc->plc3.holes[h][f].end(), hole_facet_nodes);
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_clear(polytope_plc_t* plc)
{
  if (plc->dimension == 2)
    plc->plc2.clear();
  else
    plc->plc3.clear();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
bool polytope_plc_empty(polytope_plc_t* plc)
{
  if (plc->dimension == 2)
    return plc->plc2.empty();
  else
    return plc->plc3.empty();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
bool polytope_plc_valid(polytope_plc_t* plc)
{
  if (plc->dimension == 2)
    return plc->plc2.valid();
  else
    return plc->plc3.valid();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_fprintf(polytope_plc_t* plc, FILE* stream)
{
  // FIXME
}
//------------------------------------------------------------------------

}

