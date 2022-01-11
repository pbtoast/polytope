from PYB11Generator import *
from TessellatorCommonMethods import *

@PYB11template("int Dimension", "RealType")
class Tessellator:
    """An abstract base class for objects that generate 
Voronoi and Voronoi-like tessellations for sets of points and/or 
geometries."""

    PYB11typedefs = """
  typedef typename DimensionTraits<%(Dimension)s, %(RealType)s>::QuantizedTessellation QuantizedTessellation;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"
        return

    # We provide this alias for tessellating in a box since the native polytope
    # representation of the (lower, upper) coordinates are just C-array pointers,
    # which are cumbersome for python
    @PYB11const
    @PYB11implementation("""[](Tessellator<%(Dimension)s, %(RealType)s>& self,
                               const std::vector<%(RealType)s>& points,
                               py::tuple low,
                               py::tuple high,
                               Tessellation<%(Dimension)s, %(RealType)s>& mesh) {
                                   std::vector<double> clow, chigh;
                                   for (auto i = 0; i < %(Dimension)s; ++i) {
                                     clow.push_back(low[i].cast<double>());
                                     chigh.push_back(high[i].cast<double>());
                                   }
                                   self.tessellate(points, &clow.front(), &chigh.front(), mesh);
                               }""")
    def tessellateBox(self,
                      points = "const std::vector<%(RealType)s>&",
                      low = "py::tuple",
                      high = "py::tuple",
                      mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        """Generate a Voronoi tessellation for the given set of generator points
with a bounding box specified by \a low and \a high. Here, low[i]
contains the ith coordinate for the "lower-left-near" corner of the 
bounding box in 2D or 3D, and high[i] contains the corresponding 
opposite corner. The coordinates of these points are stored in 
point-major order and the 0th component of the ith point appears in 
points[Dimension*i].
\\\\param points A (Dimension*numPoints) array containing point coordinates.
\\\\param low The coordinates of the "lower-left-near" bounding box corner.
\\\\param high The coordinates of the "upper-right-far" bounding box corner.
\\\\param mesh This will store the resulting tessellation."""
        return "void"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def tessellate(self,
                   points = "const std::vector<%(RealType)s>&",
                   mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        """Generate a Voronoi tessellation for the given set of generator points.
The coordinates of these points are stored in point-major order and 
the 0th component of the ith point appears in points[Dimension*i].
\\\\param points A (Dimension*numPoints) array containing point coordinates.
\\\\param mesh This will store the resulting tessellation."""
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("tessellate")
    def tessellate1(self,
                    points = "const std::vector<%(RealType)s>&",
                    low = "%(RealType)s*",
                    high = "%(RealType)s*",
                    mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        """Generate a Voronoi tessellation for the given set of generator points
with a bounding box specified by \a low and \a high. Here, low[i]
contains the ith coordinate for the "lower-left-near" corner of the 
bounding box in 2D or 3D, and high[i] contains the corresponding 
opposite corner. The coordinates of these points are stored in 
point-major order and the 0th component of the ith point appears in 
points[Dimension*i].
\\\\param points A (Dimension*numPoints) array containing point coordinates.
\\\\param low The coordinates of the "lower-left-near" bounding box corner.
\\\\param high The coordinates of the "upper-right-far" bounding box corner.
\\\\param mesh This will store the resulting tessellation."""
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("tessellate")
    def tessellate2(self,
                    points = "const std::vector<%(RealType)s>&",
                    PLCpoints = "const std::vector<%(RealType)s>&",
                    geometry = "const PLC<%(Dimension)s, %(RealType)s>&",
                    mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        """Generate a Voronoi-like tessellation for the given set of generator 
points and a description of the geometry in which they exist.
The coordinates of these points are stored in point-major order and 
the 0th component of the ith point appears in points[Dimension*i].
This default implementation issues an error explaining that the 
Tessellator does not support PLCs.
\\\\param points A (Dimension*numPoints) array containing point coordinates.
\\\\param PLCpoints A (Dimension*n) array containing point coordinates for the PLC.
\\\\param geometry A description of the geometry in Piecewise Linear Complex form.
\\\\param mesh This will store the resulting tessellation.
"""
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("tessellate")
    def tessellate3(self,
                    points = "const std::vector<%(RealType)s>&",
                    geometry = "const ReducedPLC<%(Dimension)s, %(RealType)s>&",
                    mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        """Generate a Voronoi-like tessellation for the given set of generator 
points and a description of the geometry in which they exist.
The geometry description uses the ReducedPLC to combine vertex
coordinates and facet topology into a single struct out of convenience.
\\\\param points A (Dimension*numPoints) array containing point coordinates.
\\\\param geometry A description of the geometry in Reduced Piecewise Linear Complex form.
\\\\param mesh This will store the resulting tessellation."""
        return "void"

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

    # @PYB11protected
    # @PYB11const
    # def computeNormalizedPoints(self,
    #                             points = "const std::vector<%(RealType)s>&",
    #                             PLCpoints = "const std::vector<%(RealType)s>&",
    #                             computeBounds = "const bool",
    #                             low = "%(RealType)s*",
    #                             high = "%(RealType)s*"):
    #     "Return a normalized set of coordinates, also returning the bounding low/high points."
    #     return "std::vector<%(RealType)s>"

#-------------------------------------------------------------------------------
# Inject the common methods
PYB11inject(TessellatorCommonMethods, Tessellator, pure_virtual=True)

#-------------------------------------------------------------------------------
# Template instantiations
Tessellator2d = PYB11TemplateClass(Tessellator, template_parameters=("2", "double"))
Tessellator3d = PYB11TemplateClass(Tessellator, template_parameters=("3", "double"))
