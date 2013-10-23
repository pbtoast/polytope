#include "polytope.hh"
#include "polytope_c.h"

using namespace polytope;

namespace 
{

static polytope_error_handler err_handler = NULL;

void polytope_c_error_handler(const std::string& message,
                              int status)
{
  err_handler(message.c_str(), status);
}

}

extern "C" 
{

void polytope_set_error_handler(polytope_error_handler handler)
{
  POLY_ASSERT(handler != NULL);
  err_handler = handler;
  setErrorHandler(polytope_c_error_handler);
}

}

