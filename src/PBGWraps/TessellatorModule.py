from pybindgen import *

import os,sys
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
        self.objs = []
        
        # --- Base Tessellators
        self.Tessellator2d = addObject(polytope, 
                                       "Tessellator2d", 
                                       allow_subclassing=True)
        self.Tessellator3d = addObject(polytope, 
                                       "Tessellator3d", 
                                       allow_subclassing=True)
        self.objs.append(self.Tessellator2d)
        self.objs.append(self.Tessellator3d)

        # --- The Triangle Tessellator
        if (mod.have_triangle):
            self.TriangleTessellator2d = addObject(polytope, 
                                                   "TriangleTessellator2d", 
                                                   parent=self.Tessellator2d)
            self.objs.append(self.TriangleTessellator2d)

        # --- The Boost.Voronoi Tessellator
        if (mod.have_boost_voronoi):
            self.BoostTessellator2d = addObject(polytope, 
                                                "BoostTessellator2d", 
                                                parent=self.Tessellator2d)
            self.objs.append(self.BoostTessellator2d)

        # --- The Tetgen Tessellator
        if (mod.have_tetgen):
            self.TetgenTessellator3d = addObject(polytope, 
                                                 "TetgenTessellator3d", 
                                                 parent=self.Tessellator3d)
            self.objs.append(self.TetgenTessellator3d)
        
        # --- The Distributed and Serial-Distributed Tessellators
        if (mod.have_mpi):
            self.DistributedTessellator2d = addObject(polytope,
                                                      "DistributedTessellator2d",
                                                      parent=self.Tessellator2d,
                                                      allow_subclassing=True)
            self.DistributedTessellator3d = addObject(polytope,
                                                      "DistributedTessellator3d",
                                                      parent=self.Tessellator3d,
                                                      allow_subclassing=True)
            self.SerialDistributedTessellator2d = addObject(polytope,
                                                            "SerialDistributedTessellator2d",
                                                            parent=self.DistributedTessellator2d)
            self.SerialDistributedTessellator3d = addObject(polytope,
                                                            "SerialDistributedTessellator3d",
                                                            parent=self.DistributedTessellator3d)
            


        # self.VoroTessellator2d = addObject(polytope, "VoroTessellator2d", parent=self.Tessellator2d)
        # self.VoroTessellator3d = addObject(polytope, "VoroTessellator3d", parent=self.Tessellator3d)
        
        return
    
    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
                       
        self.generateTessellatorBindings(self.Tessellator2d, 2)
        self.generateTessellatorBindings(self.Tessellator3d, 3)

        if (mod.have_triangle):
            self.generateTriangleTessellatorBindings(self.TriangleTessellator2d, 2)

        if (mod.have_boost_voronoi):
            self.generateBoostTessellatorBindings(self.BoostTessellator2d, 2)

        if (mod.have_tetgen):
            self.generateTetgenTessellatorBindings(self.TetgenTessellator3d, 3)
        
        if (mod.have_mpi):
            self.generateDistributedTessellatorBindings(self.DistributedTessellator2d, 2)
            self.generateDistributedTessellatorBindings(self.DistributedTessellator3d, 3)
            self.generateSerialDistributedTessellatorBindings(self.SerialDistributedTessellator2d, 2)
            self.generateSerialDistributedTessellatorBindings(self.SerialDistributedTessellator3d, 3)

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
                                          param("double *", "low"),
                                          param("double *", "high"),
                                          refparam(Tessellation, "mesh")],
                     is_virtual=True, is_const=True)
        x.add_method("tessellate", None, [constrefparam("vector_of_double", "points"),
                                          constrefparam("vector_of_double", "PLCpoints"),
                                          constrefparam(PLC, "geometry"),
                                          refparam(Tessellation, "mesh")],
                     is_virtual=True, is_const=True)

        x.add_method("tessellateDegenerate", "vector_of_unsigned", 
                     [constrefparam("vector_of_double", "points"),
                      param("double", "tol"),
                      refparam(Tessellation, "mesh")],
                     is_virtual=True, is_const=True)
        x.add_method("tessellateDegenerate", "vector_of_unsigned", 
                     [constrefparam("vector_of_double", "points"),
                      param("double *", "low"),
                      param("double *", "high"),
                      param("double", "tol"),
                      refparam(Tessellation, "mesh")],
                     is_virtual=True, is_const=True)
        x.add_method("tessellateDegenerate", "vector_of_unsigned",
                     [constrefparam("vector_of_double", "points"),
                      constrefparam("vector_of_double", "PLCpoints"),
                      constrefparam(PLC, "geometry"),
                      param("double", "tol"),
                      refparam(Tessellation, "mesh")],
                     is_virtual=True, is_const=True)

        x.add_method("handlesPLCs", retval("bool"), [], is_pure_virtual=True, is_const=True)
        x.add_method("name", retval("std::string"), [], is_pure_virtual=True, is_const=True)
        x.add_method("degeneracy", "double", [], is_pure_virtual=True, is_const=True)

        return

    #---------------------------------------------------------------------------
    # Bindings (TriangleTessellator)
    #---------------------------------------------------------------------------
    def generateTriangleTessellatorBindings(self, x, ndim):

        # Constructors
        x.add_constructor([])

        # Methods
        x.add_method("handlesPLCs", retval("bool"), [], is_virtual=True, is_const=True)
        x.add_method("name", retval("std::string"), [], is_virtual=True, is_const=True)
        x.add_method("degeneracy", "double", [], is_virtual=True, is_const=True)
        
        return

    #---------------------------------------------------------------------------
    # Bindings (BoostTessellator)
    #---------------------------------------------------------------------------
    def generateBoostTessellatorBindings(self, x, ndim):

        # Constructors
        x.add_constructor([])

        # Methods
        x.add_method("handlesPLCs", retval("bool"), [], is_virtual=True, is_const=True)
        x.add_method("name", retval("std::string"), [], is_virtual=True, is_const=True)
        x.add_method("degeneracy", "double", [], is_virtual=True, is_const=True)
        
        return

    #---------------------------------------------------------------------------
    # Bindings (TetgenTessellator)
    #---------------------------------------------------------------------------
    def generateTetgenTessellatorBindings(self, x, ndim):

        # Constructors
        x.add_constructor([])

        # Methods
        x.add_method("handlesPLCs", retval("bool"), [], is_virtual=True, is_const=True)
        x.add_method("name", retval("std::string"), [], is_virtual=True, is_const=True)
        x.add_method("degeneracy", "double", [], is_virtual=True, is_const=True)
        
        return    
    
    #---------------------------------------------------------------------------
    # Bindings (DistributedTessellator)
    #---------------------------------------------------------------------------
    def generateDistributedTessellatorBindings(self, x, ndim):
        
        # Object names
        TessellatorPtr  = "polytope::Tessellator%id*" % ndim
        
        # Constructors
        x.add_constructor([param(TessellatorPtr, "serialTessellator", transfer_ownership=False),
                           param("bool", "assumeControl", default_value="true"),
                           param("bool", "buildCommunicationInfo", default_value="true")])
        
        # Methods
        x.add_method("handlesPLCs", retval("bool"), [], is_virtual=True, is_const=True)
        x.add_method("name", retval("std::string"), [], is_virtual=True, is_const=True)
        x.add_method("degeneracy", "double", [], is_virtual=True, is_const=True)
        
        return

    #---------------------------------------------------------------------------
    # Bindings (SerialDistributedTessellator)
    #---------------------------------------------------------------------------
    def generateSerialDistributedTessellatorBindings(self, x, ndim):
        
        # Object names
        Tessellator  = "polytope::Tessellator%id" % ndim
        
        # Constructors
        x.add_constructor([inputptrparam(Tessellator, "serialTessellator"),
                           param("bool", "assumeControl", default_value="true"),
                           param("bool", "buildCommunicationInfo", default_value="true")])
                
        return

    """
    #---------------------------------------------------------------------------
    # Bindings (VoroTessellator)
    #---------------------------------------------------------------------------
    def generateTriangleTessellatorBindings(self, x, ndim):
        return
    """
    
    
