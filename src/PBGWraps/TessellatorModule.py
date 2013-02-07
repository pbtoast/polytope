from pybindgen import *

import sys
sys.path.append(".")
from PBGutils import *

#-------------------------------------------------------------------------------
# Class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Tessellator:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):
        
        # Includes
        mod.add_include('"PolytopeTypes.hh"')

        # Namespaces
        polytope = mod.add_cpp_namespace("polytope")

        # Expose types
        self.Tessellator2d = addObject(polytope, "Tessellator2d", allow_subclassing=True)
        self.Tessellator3d = addObject(polytope, "Tessellator3d", allow_subclassing=True)

        self.TriangleTessellator2d = addObject(polytope, "TriangleTessellator2d", parent=self.Tessellator2d)
        
        self.TetgenTessellator3d = addObject(polytope, "TetgenTessellator3d", parent=self.Tessellator3d)
        
        # self.VoroTessellator2d = addObject(polytope, "VoroTessellator2d", parent=self.Tessellator2d)
        # self.VoroTessellator3d = addObject(polytope, "VoroTessellator3d", parent=self.Tessellator3d)

        return
    
    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
                       
        self.generateTessellatorBindings(self.Tessellator2d, 2)
        self.generateTessellatorBindings(self.Tessellator3d, 3)

        self.generateTriangleTessellatorBindings(self.TriangleTessellator2d, 2)

        self.generateTetgenTessellatorBindings(self.TetgenTessellator3d, 3)
        
        # self.generateVoroTessellatorBindings(self.VoroTessellator2d, 2)
        # self.generateVoroTessellatorBindings(self.VoroTessellator3d, 3)

        return

    #---------------------------------------------------------------------------
    # Bindings (Tessellator)
    #---------------------------------------------------------------------------
    def generateTessellatorBindings(self, x, ndim):
        
        # Object names
        PLC          = "polytope::PLC%id" % ndim
        Tessellation = "polytope::Tessellation%id" % ndim
        
        # Constructors
        x.add_constructor([])

        # Methods
        x.add_method("tessellate", None, [constrefparam("vector_of_double", "points"),
                                          refparam(Tessellation, "mesh")],
                     is_virtual=True, is_const=True)
        x.add_method("tessellate", None, [constrefparam("vector_of_double", "points"),
                                          param("double *", "low"),#, transfer_ownership=True),
                                          param("double *", "high"),#, transfer_ownership=True),
                                          refparam(Tessellation, "mesh")],
                     is_virtual=True, is_const=True)
        x.add_method("tessellate", None, [constrefparam("vector_of_double", "points"),
                                          constrefparam("vector_of_double", "PLCpoints"),
                                          constrefparam(PLC, "geometry"),
                                          refparam(Tessellation, "mesh")],
                     is_virtual=True, is_const=True)
        x.add_method("handlesPLCs", retval('bool'), [], is_pure_virtual=True, is_const=True)
        
        # Attributes
        
        
        return

    #---------------------------------------------------------------------------
    # Bindings (TriangleTessellator)
    #---------------------------------------------------------------------------
    def generateTriangleTessellatorBindings(self, x, ndim):

        # Constructors
        x.add_constructor([])

        # Methods
        x.add_method("handlesPLCs", retval("bool"), [], is_virtual=True, is_const=True)
        
        return

    #---------------------------------------------------------------------------
    # Bindings (TetgenTessellator)
    #---------------------------------------------------------------------------
    def generateTetgenTessellatorBindings(self, x, ndim):

        # Constructors
        x.add_constructor([param("const bool", "directComputation")])

        # Methods
        x.add_method("handlesPLCs", retval("bool"), [], is_virtual=True, is_const=True)
        x.add_method("directComputation", retval("bool"), [], is_const=True)
        x.add_method("directComputation", None, [param("const bool", "x")])

        return
    
    
    """
    #---------------------------------------------------------------------------
    # Bindings (VoroTessellator)
    #---------------------------------------------------------------------------
    def generateTriangleTessellatorBindings(self, x, ndim):
        return
    """
    
    
