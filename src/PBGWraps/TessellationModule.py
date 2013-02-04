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
        x.add_method("computeNodeCells", retval("std::vector<std::set<unsigned> >"), [])
        x.add_method("computeCellToNodes", retval("std::vector<std::set<unsigned> >"), [])
        
        # Attributes
        x.add_custom_instance_attribute("nodes", retval("std::vector<double>*", caller_owns_return=False), getter="polytope::getNodes", setter="polytope::setNodes", getter_template_parameters=[str(ndim),"double"], setter_template_parameters=[str(ndim),"double"])
        x.add_instance_attribute("cells", "std::vector<std::vector<int> >")
        x.add_instance_attribute("faces", "std::vector<std::vector<unsigned> >")
        x.add_instance_attribute("faceCells", "std::vector<std::vector<int> >")
        x.add_instance_attribute("convexHull", PLC)
        x.add_instance_attribute("neighborDomains","std::vector<unsigned>")
        x.add_instance_attribute("sharedNodes","std::vector<std::vector<unsigned> >")
        x.add_instance_attribute("sharedFaces","std::vector<std::vector<unsigned> >")
        
        return
