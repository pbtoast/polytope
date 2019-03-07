from PYB11Generator import *
from TessellatorCommonMethods import *

@PYB11template("int Dimension", "RealType")
class Tessellator:
    """An abstract base class for objects that generate 
Voronoi and Voronoi-like tessellations for sets of points and/or 
geometries."""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"
        return

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def tessellateDegenerate(self,
                             points = "const std::vector<%(RealType)s>&",
                             tol = "const %(RealType)s",
                             mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        return "std::vector<unsigned>"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("tessellateDegenerate")
    def tessellateDegenerate1(self,
                              points = "const std::vector<%(RealType)s>&",
                              low = "%(RealType)s*",
                              high = "%(RealType)s*",
                              tol = "const %(RealType)s",
                              mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        return "std::vector<unsigned>"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("tessellateDegenerate")
    def tessellateDegenerate2(self,
                              points = "const std::vector<%(RealType)s>&",
                              PLCpoints = "const std::vector<%(RealType)s>&",
                              geometry = "const PLC<%(Dimension)s, %(RealType)s>&",
                              tol = "const %(RealType)s",
                              mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        return "std::vector<unsigned>"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("tessellateDegenerate")
    def tessellateDegenerate3(self,
                              points = "const std::vector<%(RealType)s>&",
                              geometry = "const ReducedPLC<%(Dimension)s, %(RealType)s>&",
                              tol = "const %(RealType)s",
                              mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        return "std::vector<unsigned>"

    @PYB11virtual
    @PYB11const
    def tessellateNormalized(self,
                             points = "const std::vector<%(RealType)s>&",
                             mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("tessellateNormalized")
    def tessellateNormalized1(self,
                              points = "const std::vector<%(RealType)s>&",
                              low = "%(RealType)s*",
                              high = "%(RealType)s*",
                              mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("tessellateNormalized")
    def tessellateNormalized2(self,
                              points = "const std::vector<%(RealType)s>&",
                              PLCpoints = "const std::vector<%(RealType)s>&",
                              geometry = "const PLC<%(Dimension)s, %(RealType)s>&",
                              mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("tessellateNormalized")
    def tessellateNormalized3(self,
                              points = "const std::vector<%(RealType)s>&",
                              geometry = "const ReducedPLC<%(Dimension)s, %(RealType)s>&",
                              mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        return "void"

    @PYB11pure_virtual
    @PYB11const
    def handlesPLCs(self):
        """Override this method to return true if this Tessellator supports 
the description of a domain boundary using a PLC (as in the second 
tessellate method, above), and false if it does not. Some algorithms 
for tessellation do not naturally accommodate an explicit boundary 
description, and Tessellators using these algorithms should override 
this method to return false. A stub method for PLC-enabled
tessellation is provided for convenience.
This query mechanism prevents us from descending into the taxonomic 
hell associated with elaborate inheritance hierarchies."""
        return "bool"

    @PYB11pure_virtual
    @PYB11const
    def name(self):
        "A unique name string per tessellation instance."
        return "std::string"

    @PYB11pure_virtual
    @PYB11const
    def degeneracy(self):
        """Returns the accuracy to which this tessellator can distinguish coordinates.
Should be returned appropriately for normalized coordinates, i.e., if all
coordinates are in the range xi \\\\in [0,1], what is the minimum allowed 
delta in x."""
        return "%(RealType)s"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("degeneracy")
    def setdegeneracy(self,
                      val = "const %(RealType)s"):
        return "void"

    #...........................................................................
    # Protected methods
#     @PYB11protected
#     @PYB11const
#     def boundingBox(self,
#                     points = "std::vector<%(RealType)s>&"):
#         """This helper method creates a piecewise linear complex (PLC) 
# representing the bounding box containing the given points and 
# adds the corners of the bounding box to \a points."""
#         return "PLC<%(Dimension)s, %(RealType)s>"

    @PYB11protected
    @PYB11const
    def computeNormalizedPoints(self,
                                points = "const std::vector<%(RealType)s>&",
                                PLCpoints = "const std::vector<%(RealType)s>&",
                                computeBounds = "const bool",
                                low = "%(RealType)s*",
                                high = "%(RealType)s*"):
        "Return a normalized set of coordinates, also returning the bounding low/high points."
        return "std::vector<%(RealType)s>"

#-------------------------------------------------------------------------------
# Inject the common methods
PYB11inject(TessellatorCommonMethods, Tessellator, virtual=True)

#-------------------------------------------------------------------------------
# Template instantiations
Tessellator2d = PYB11TemplateClass(Tessellator, template_parameters=("2", "double"))
Tessellator3d = PYB11TemplateClass(Tessellator, template_parameters=("3", "double"))
