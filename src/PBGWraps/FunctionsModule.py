from pybindgen import *

import sys
sys.path.append(".")
from PBGutils import *

#-------------------------------------------------------------------------------
# Class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Functions:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):
        
        # Includes
        mod.add_include('"PolytopeTypes.hh"')

        # Namespaces
        polytope = mod.add_cpp_namespace("polytope")

        self.objs = []

        return
    
    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        polytope = mod.add_cpp_namespace("polytope")
                       
        polytope.add_function("writePLCtoOFF",
                              None,
                              [constrefparam("polytope::PLC3d", "plc"),
                               constrefparam("vector_of_double", "coords"),
                               param("std::string", "filename")],
                              template_parameters = ["double"],
                              custom_name = "writePLCtoOFF",
                              docstring = "writePLCtoOFF -- output a PLC to an OFF file for use with GeomView.")


        for dim in (2, 3):
            polytope.add_function("writeTessellation", 
                                  None,
                                  [constrefparam("polytope::Tessellation%id" % dim, "mesh"),
                                   param("std::string", "filePrefix"),
                                   param("PyObject*", "nodeFields", transfer_ownership=False),
                                   param("PyObject*", "edgeFields", transfer_ownership=False),
                                   param("PyObject*", "faceFields", transfer_ownership=False),
                                   param("PyObject*", "cellFields", transfer_ownership=False),
                                   param("int", "cycle", default_value = "0"),
                                   param("double", "time", default_value = "0.0")],
                                  template_parameters = [str(dim), "double"],
                                  custom_name = ("writeTessellation%id" % dim),
                                  docstring = "writeTessellation -- output a silo file with optional fields")
# arguments:
#   mesh       : the tessellation to write
#   filePrefix : prefix for the output file name
#   nodeFields : (optional, default=None) dictionary of (name -> vals) node centered values to put in file
#   edgeFields : (optional, default=None) dictionary of (name -> vals) edge centered values to put in file
#   faceFields : (optional, default=None) dictionary of (name -> vals) face centered values to put in file
#   cellFields : (optional, default=None) dictionary of (name -> vals) cell centered values to put in file
#   cycle      : (optional, default=0) cycle number
#   time       : (optional, default=0.0) time
# """,

            polytope.add_function("writeTessellation", 
                                  None,
                                  [constrefparam("polytope::Tessellation%id" % dim, "mesh"),
                                   param("std::string", "filePrefix"),
                                   param("std::string", "directory"),
                                   param("PyObject*", "nodeFields", transfer_ownership=False),
                                   param("PyObject*", "edgeFields", transfer_ownership=False),
                                   param("PyObject*", "faceFields", transfer_ownership=False),
                                   param("PyObject*", "cellFields", transfer_ownership=False),
                                   param("int", "cycle", default_value = "0"),
                                   param("double", "time", default_value = "0.0")],
                                  template_parameters = [str(dim), "double"],
                                  custom_name = ("writeTessellation%id" % dim),
                                  docstring = "writeTessellation -- output a silo file with optional fields")
# arguments:
#   mesh       : the tessellation to write
#   filePrefix : prefix for the output file name
#   directory  : the base directory for the files
#   nodeFields : (optional, default=None) dictionary of (name -> vals) node centered values to put in file
#   edgeFields : (optional, default=None) dictionary of (name -> vals) edge centered values to put in file
#   faceFields : (optional, default=None) dictionary of (name -> vals) face centered values to put in file
#   cellFields : (optional, default=None) dictionary of (name -> vals) cell centered values to put in file
#   cycle      : (optional, default=0) cycle number
#   time       : (optional, default=0.0) time
# """,


            polytope.add_function("constructConvexHull%id" % dim,
                                  retval("polytope::PLC%id*" % dim, caller_owns_return=True),
                                  [constrefparam("vector_of_double", "points"),
                                   constrefparam("vector_of_double", "points"),
                                   param("double", "dx")],
                                  template_parameters = ["double"],
                                  custom_name = ("constructConvexHull%id" % dim),
                                  docstring = "constructConvexHull -- output the convex hull of a collection of points as a PLC")


        return
