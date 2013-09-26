#include "polytope_c.h"
#include "polytope.hh"
#include "BoostTessellator.hh"
#include "TriangleTessellator.hh"
#include "TetgenTessellator.hh"

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
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellator_tessellate_in_box(polytope_tessellator_t* tessellator,
                                            polytope_real_t* points, int num_points,
                                            polytope_real_t* low,
                                            polytope_real_t* high,
                                            polytope_tessellation_t* mesh)
{
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellator_tessellate_in_plc(polytope_tessellator_t* tessellator,
                                            polytope_real_t* points, int num_points,
                                            polytope_real_t* plc_points, int num_plc_points,
                                            polytope_plc_t* piecewise_linear_complex,
                                            polytope_tessellation_t* mesh)
{
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
#if HAVE_BOOST
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
#if HAVE_TRIANGLE
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
#if HAVE_TETGEN
polytope_tessellator_t* tetgen_tessellator_new()
{
  polytope_tessellator_t* t = (polytope_tessellator_t*)malloc(sizeof(polytope_tessellator_t));
  t->tess2 = NULL;
  t->tess3 = new TetgenTessellator();
  return t;
}
#endif
//------------------------------------------------------------------------

//------------------------------------------------------------------------
polytope_tessellator_t* voroplusplus_tessellator_new(int dimension)
{
  polytope_tessellator_t* t = (polytope_tessellator_t*)malloc(sizeof(polytope_tessellator_t));
  if (dimension == 2)
  {
    t->tess2 = new VoroPP_2d<polytope_real_t>();
    t->tess3 = NULL;
  }
  else
  {
    t->tess2 = NULL;
    t->tess3 = new VoroPP_3d<polytope_real_t>();
  }
  return t;
}
//------------------------------------------------------------------------

}

