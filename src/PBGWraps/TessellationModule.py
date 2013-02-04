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

        return
    
    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
                       
        self.generateTessellationBindings(self.Tessellation2d, 2)
        self.generateTessellationBindings(self.Tessellation3d, 3)

        return

    #---------------------------------------------------------------------------
    # Bindings (Tessellation)
    #---------------------------------------------------------------------------
    def generateTessellationBindings(self, x, ndim):
        
        # Object names
        PLC = "polytope::PLC%id" % ndim
        
        # Constructors
        x.add_constructor([])
        
        # Methods
        x.add_method("clear", None, [])
        x.add_method("empty", retval("bool"), [])
        x.add_method("computeNodeCells", retval("vector_of_set_of_unsigned"), [])
        x.add_method("computeCellToNodes", retval("vector_of_set_of_unsigned"), [])
        
        # Attributes
        x.add_custom_instance_attribute("nodes", retval("vector_of_double*", caller_owns_return=False), getter="polytope::getNodes", setter="polytope::setNodes", getter_template_parameters=[str(ndim),"double"], setter_template_parameters=[str(ndim),"double"])
        x.add_instance_attribute("cells", "vector_of_vector_of_int")
        x.add_instance_attribute("faces", "vector_of_vector_of_unsigned")
        x.add_instance_attribute("faceCells", "vector_of_vector_of_int")
        x.add_instance_attribute("convexHull", PLC)
        x.add_instance_attribute("neighborDomains","vector_of_unsigned")
        x.add_instance_attribute("sharedNodes","vector_of_vector_of_unsigned")
        x.add_instance_attribute("sharedFaces","vector_of_vector_of_unsigned")
        
        return
