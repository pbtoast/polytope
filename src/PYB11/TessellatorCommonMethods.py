from PYB11Generator import *

@PYB11ignore
class TessellatorCommonMethods:

    #...........................................................................
    # Virtual methods
    @PYB11const
    def tessellate(self,
                   points = "const std::vector<%(RealType)s>&",
                   mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        """Generate a Voronoi tessellation for the given set of generator points.
The coordinates of these points are stored in point-major order and 
the 0th component of the ith point appears in points[Dimension*i].
\param points A (Dimension*numPoints) array containing point coordinates.
\param mesh This will store the resulting tessellation."""
        return "void"

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
\param points A (Dimension*numPoints) array containing point coordinates.
\param low The coordinates of the "lower-left-near" bounding box corner.
\param high The coordinates of the "upper-right-far" bounding box corner.
\param mesh This will store the resulting tessellation."""
        return "void"

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
\param points A (Dimension*numPoints) array containing point coordinates.
\param PLCpoints A (Dimension*n) array containing point coordinates for the PLC.
\param geometry A description of the geometry in Piecewise Linear Complex form.
\param mesh This will store the resulting tessellation.
"""
        return "void"

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
\param points A (Dimension*numPoints) array containing point coordinates.
\param geometry A description of the geometry in Reduced Piecewise Linear Complex form.
\param mesh This will store the resulting tessellation."""
        return "void"

