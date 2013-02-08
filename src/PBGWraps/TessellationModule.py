from pybindgen import *

import sys
sys.path.append(".")
from PBGutils import *

#-------------------------------------------------------------------------------
# Class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Tessellation:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):
        
        # Includes
        mod.add_include('"PolytopeTypes.hh"')

        # Namespaces
        polytope = mod.add_cpp_namespace("polytope")

        # Expose types
        self.Tessellation2d = addObject(polytope, "Tessellation2d")
        self.Tessellation3d = addObject(polytope, "Tessellation3d")

        self.objs = [self.Tessellation2d, self.Tessellation3d]

        return
    
    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
                       
        for (obj, dim) in ((self.Tessellation2d, 2),
                           (self.Tessellation3d, 3)):
            self.generateTessellationBindings(obj, dim)

        return

    #---------------------------------------------------------------------------
    # Bindings (Tessellation)
    #---------------------------------------------------------------------------
    def generateTessellationBindings(self, x, ndim):
        
        # Object names
        PLC    = "polytope::PLC%id"  % ndim
        PLCptr = "polytope::PLC%id*" % ndim
        
        # Constructors
        x.add_constructor([])
        
        # Methods
        x.add_method("clear", None, [])
        x.add_method("empty", retval("bool"), [])
        x.add_method("computeNodeCells", retval("vector_of_set_of_unsigned"), [])
        x.add_method("computeCellToNodes", retval("vector_of_set_of_unsigned"), [])
        
        # Attributes
        attributes = ["nodes", "cells", "faces", "infNodes", "faceCells",
                      "convexHull", "neighborDomains", "sharedNodes", "sharedFaces"]

        returnvals = ["vector_of_double*",
                      "vector_of_vector_of_int*", 
                      "vector_of_vector_of_unsigned*", 
                      "vector_of_unsigned*",
                      "vector_of_vector_of_int*", 
                      PLCptr,
                      "vector_of_unsigned*", 
                      "vector_of_vector_of_unsigned*", 
                      "vector_of_vector_of_unsigned*"]
        
        for i,(att,ret) in enumerate(zip(attributes,returnvals)):
            x.add_custom_instance_attribute(att, retval(ret, reference_existing_object=True),
                                            getter="polytope::get"+att,
                                            setter="polytope::set"+att,
                                            getter_template_parameters=[str(ndim),"double"], 
                                            setter_template_parameters=[str(ndim),"double"])
                
        return
