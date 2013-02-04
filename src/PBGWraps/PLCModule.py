from pybindgen import *

import sys
sys.path.append(".")
from PBGutils import *

#-------------------------------------------------------------------------------
# Class to handle wrapping this module.
#-------------------------------------------------------------------------------
class PLC:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):
        
        # Includes
        mod.add_include('"PolytopeTypes.hh"')

        # Namespaces
        polytope = mod.add_cpp_namespace("polytope")

        # Expose types
        self.PLC2d = addObject(polytope, "PLC2d")
        self.PLC3d = addObject(polytope, "PLC3d")

        return
    
    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
                       
        self.generatePLCBindings(self.PLC2d, 2)
        self.generatePLCBindings(self.PLC3d, 3)

        return

    #---------------------------------------------------------------------------
    # Bindings (PLC)
    #---------------------------------------------------------------------------
    def generatePLCBindings(self, x, ndim):
        
        # Constructors
        x.add_constructor([])
        
        # Methods
        x.add_method("clear", None, [])
        x.add_method("empty", retval("bool"), [], is_const=True)
        x.add_method("valid", retval("bool"), [], is_const=True)
        
        # Attributes
        x.add_instance_attribute("facets", "std::vector<std::vector<int> >")
        x.add_instance_attribute("holes", "std::vector<std::vector<std::vector<int> > >")
                
        return
