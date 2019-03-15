from PYB11Generator import *
from Tessellator import *
from TessellatorCommonMethods import *

@PYB11template("RealType")                                 # Override base class template parameters
@PYB11template_dict({"Dimension" : "2"})
class TriangleTessellator(Tessellator):
    """An implemenation of the Tessellator interface that uses the Triangle
library by Jonathan Shewchuk."""

    PYB11typedefs = """
typedef typename DimensionTraits<%(Dimension)s, %(RealType)s>::QuantizedTessellation QuantizedTessellation;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"
        return

    #...........................................................................
    # Virtual methods
    @PYB11virtual
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

    @PYB11virtual
    @PYB11const
    def name(self):
        "A unique name string per tessellation instance."
        return "std::string"

    @PYB11virtual
    @PYB11const
    def degeneracy(self):
        """Returns the accuracy to which this tessellator can distinguish coordinates.
Should be returned appropriately for normalized coordinates, i.e., if all
coordinates are in the range xi \\\\in [0,1], what is the minimum allowed 
delta in x."""
        return "%(RealType)s"

#-------------------------------------------------------------------------------
# Inject the common methods
PYB11inject(TessellatorCommonMethods, TriangleTessellator)

#-------------------------------------------------------------------------------
# Template instantiations
TriangleTessellator2d = PYB11TemplateClass(TriangleTessellator, template_parameters="double")
