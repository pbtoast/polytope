#include "polytope_plc.h"
#include "polytope.hh"

using namespace std;
using namespace polytope;

#define MY_PLC(plc_struct) ((plc_struct->plc2 != NULL) ? plc_struct->plc2 : plc_struct->plc3)

extern "C"
{

struct polytope_plc_t
{
  PLC<2, polytope_real_t>* plc2;
  PLC<3, polytope_real_t>* plc3;
}

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
  delete MY_PLC(plc);
  free(plc);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_add_facet(polytope_plc_t* plc)
{
  MY_PLC(plc)->facets.push_back(vector<int>());
  return MY_PLC(plc)->facets.size()-1;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_add_facet_node(polytope_plc_t* plc, int facet, int node)
{
  size_t f = (size_t)facet;
  POLY_ASSERT(f < MY_PLC(plc)->facets.size());
  MY_PLC(plc)->facets[f].push_back(node);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_add_hole(polytope_plc_t* plc)
{
  MY_PLC(plc)->holes.push_back(vector<int>());
  return MY_PLC(plc)->holes.size()-1;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_add_hole_facet(polytope_plc_t* plc, int hole)
{
  size_t h = (size_t)hole;
  POLY_ASSERT(h < MY_PLC(plc)->holes.size());
  MY_PLC(plc)->holes[h].push_back(vector<int>());
  return MY_PLC(plc)->holes[h].size()-1;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_add_hole_facet_node(polytope_plc_t* plc, int facet, int node)
{
  size_t h = (size_t)hole;
  POLY_ASSERT(h < MY_PLC(plc)->holes.size());
  size_t f = (size_t)facet;
  POLY_ASSERT(f < MY_PLC(plc)->holes[h].size());
  MY_PLC(plc)->holes[h][f].push_back(node);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_num_facet_nodes(polytope_plc_t* plc, int facet)
{
  size_t f = (size_t)facet;
  POLY_ASSERT(f < MY_PLC(plc)->facets.size());
  return (int)(MY_PLC(plc)->facets[f].size());
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_get_facet_nodes(polytope_plc_t* plc, int facet, int* facet_nodes)
{
  size_t f = (size_t)facet;
  POLY_ASSERT(f < MY_PLC(plc)->facets.size());
  copy(MY_PLC(plc)->facets.begin(), MY_PLC(plc)->facets.end(), facet_nodes);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_num_hole_facets(polytope_plc_t* plc, int hole)
{
  size_t h = (size_t)hole;
  POLY_ASSERT(h < MY_PLC(plc)->holes.size());
  return (int)(MY_PLC(plc)->holes[h].size());
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_plc_num_hole_facet_nodes(polytope_plc_t* plc, int hole, int hole_facet)
{
  size_t h = (size_t)hole;
  POLY_ASSERT(h < MY_PLC(plc)->holes.size());
  size_t f = (size_t)hole_facet;
  POLY_ASSERT(f < MY_PLC(plc)->holes[h].size());
  return (int)(MY_PLC(plc)->holes[h][f].size());
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_get_hole_facet_nodes(polytope_plc_t* plc, int hole, int hole_facet, int* hole_facet_nodes)
{
  size_t h = (size_t)hole;
  POLY_ASSERT(h < MY_PLC(plc)->holes.size());
  size_t f = (size_t)hole_facet;
  POLY_ASSERT(f < MY_PLC(plc)->plc2.holes[h].size());
  copy(MY_PLC(plc)->holes[h][f].begin(), MY_PLC(plc)->holes[h][f].end(), hole_facet_nodes);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_clear(polytope_plc_t* plc)
{
  MY_PLC(plc)->clear();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
bool polytope_plc_empty(polytope_plc_t* plc)
{
  return MY_PLC(plc)->empty();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
bool polytope_plc_valid(polytope_plc_t* plc)
{
  return MY_PLC(plc)->valid();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_plc_fprintf(polytope_plc_t* plc, FILE* stream)
{
  // FIXME
}
//------------------------------------------------------------------------

}

