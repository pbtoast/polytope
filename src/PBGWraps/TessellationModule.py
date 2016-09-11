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
        self.QuantizedTessellation2d = addObject(polytope, "QuantizedTessellation2d",
                                                 template_parameters = ["int", "double"],
                                                 custom_name = "QuantizedTessellation2d")
        self.QuantizedTessellation3d = addObject(polytope, "QuantizedTessellation3d",
                                                 template_parameters = ["int", "double"],
                                                 custom_name = "QuantizedTessellation2d")

        self.objs = [self.Tessellation2d, self.Tessellation3d,
                     self.QuantizedTessellation2d, self.QuantizedTessellation3d]

        return
    
    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
                       
        for (obj, dim) in ((self.Tessellation2d, 2),
                           (self.Tessellation3d, 3)):
            self.generateTessellationBindings(obj, dim)
        for (obj, dim) in ((self.QuantizedTessellation2d, 2),
                           (self.QuantizedTessellation3d, 3)):
            self.generateQuantizedTessellationBindings(obj, dim)

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
        for (att,ret) in [("nodes",           "vector_of_double*"),              
                          ("cells",           "vector_of_vector_of_int*"),       
                          ("faces",           "vector_of_vector_of_unsigned*"),  
                          ("boundaryNodes",   "vector_of_unsigned*"),            
                          ("boundaryFaces",   "vector_of_unsigned*"),            
                          ("faceCells",       "vector_of_vector_of_int*"),       
                          ("convexHull",      PLCptr),                           
                          ("neighborDomains", "vector_of_unsigned*"),            
                          ("sharedNodes",     "vector_of_vector_of_unsigned*"),  
                          ("sharedFaces",     "vector_of_vector_of_unsigned*")]:
            x.add_custom_instance_attribute(att, retval(ret, reference_existing_object=True),
                                            getter="polytope::get"+att,
                                            setter="polytope::set"+att,
                                            getter_template_parameters=[str(ndim),"double"], 
                                            setter_template_parameters=[str(ndim),"double"])
                
        return

    #---------------------------------------------------------------------------
    # Bindings (QuantizedTessellation)
    #---------------------------------------------------------------------------
    def generateQuantizedTessellationBindings(self, x, ndim):
        
        # Object names
        
        # Constructors
        x.add_constructor([param("vector_of_double", "points"),
                           param("vector_of_double", "boundaryPoints")])
        
        # Methods
        
        # Attributes
                
        return
