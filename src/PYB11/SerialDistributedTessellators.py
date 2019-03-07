from PYB11Generator import *
from DistributedTessellators import *
from TessellatorCommonMethods import *

@PYB11template("int Dimension", "RealType")
class SerialDistributedTessellator(DistributedTessellator):
    """SerialDistributedTessellator

Provides a parallel tessellation.

This class assumes the user provides:
1.  A serial Tessellator.
2.  The generators in parallel, distributed in any manner the user likes
    so long as the generators are not degenerate: i.e., don't repeast the
    same generator on different domains.

This is a dumb but reliable version of DistributedTessellator.  This version
exchanges all generators to everyone, builds the full tessellation on each
processor, and then just keeps the local portion needed.  Obviously not 
scalable, but useful as a check for the real DistributedTessellator.
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               serialTessellator = "Tessellator<%(Dimension)s, %(RealType)s>*",
               assumeControl = ("bool", "true"),
               buildCommunicationInfo = ("bool", "false")):
        return

#-------------------------------------------------------------------------------
# Inject the common methods
PYB11inject(TessellatorCommonMethods, SerialDistributedTessellator)

#-------------------------------------------------------------------------------
# Template instantiations
SerialDistributedTessellator2d = PYB11TemplateClass(SerialDistributedTessellator, template_parameters=("2", "double"))
SerialDistributedTessellator3d = PYB11TemplateClass(SerialDistributedTessellator, template_parameters=("3", "double"))
