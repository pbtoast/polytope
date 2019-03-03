from PYB11Generator import *

@PYB11template("Dimension", "RealType")
class PLC:
    """class PLC - A Piecewise Linear Complex in 3D, or a Planar Straight Line 
Graph (PSLG) in 2D."""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Methods
    def clear(self):
        "Clears facets and holes to empty the PLC"
        return "void"

    @PYB11const
    def empty(self):
        "Returns true if this PLC is empty, false otherwise."
        return "bool"

    @PYB11const
    def valid(self):
        """Returns true if this PLC is valid (at first glance), false if it 
is obviously invalid. This is not a rigorous check!"""
        return "bool"

    @PYB11implementation("[](const PLC<%(Dimension)s, %(RealType)s>& self) { std::stringstream ss; ss << self; return ss.str();}")
    def __str__(self):
        return "std::string"

    #...........................................................................
    # Attributes
    facets = PYB11readwrite(doc="""This two-dimensional array defines the topology of the facets of the 
piecewise linear complex in terms of connections to generating points. 
A facet has an arbitrary number of points in 3D and 2 points in 2D. 
facets[i][j] gives the index of the jth generating point of the ith 
facet.""")
    holes = PYB11readwrite(doc="""This three dimensional array defines the topology of the inner facets
or holes in the geometry.  The outer-most dimension is the number of 
holes, and the remaining are facets using the same convention as the
the "facets" member.  In other words, holes[k][i][j] is the jth
generating point of the ith facet of the kth hole.""")

#-------------------------------------------------------------------------------
# Template instantiations
PLC2d = PYB11TemplateClass(PLC, template_parameters=("2", "double"))
PLC3d = PYB11TemplateClass(PLC, template_parameters=("3", "double"))
