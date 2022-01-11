from PYB11Generator import *
from PLC import PLC

@PYB11template("Dimension", "RealType")
class ReducedPLC(PLC):
    """class ReducedPLC - A Piecewise Linear Complex in 3D, or a Planar Straight Line 
raph (PSLG) in 2D.
he reduced PLC is a PLC which
a) contains its own generating points,
b) is reduced to just the generating points that are used in the PLC."""

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                plc = "const PLC<%(Dimension)s, %(RealType)s>&",
                allpoints = "const std::vector<%(RealType)s>&"):
        "Construct from a normal PLC, copying the necessary generator coordinates to our internal data."

    #...........................................................................
    # Attributes
    points = PYB11readwrite(doc="""This array of size (Dimension*numPoints) contains components of 
the points that are used in the facets of the PLC.""")

#-------------------------------------------------------------------------------
# Template instantiations
ReducedPLC2d = PYB11TemplateClass(ReducedPLC, template_parameters=("2", "double"))
ReducedPLC3d = PYB11TemplateClass(ReducedPLC, template_parameters=("3", "double"))
