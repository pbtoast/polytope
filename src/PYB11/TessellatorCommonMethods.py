from PYB11Generator import *

@PYB11ignore
class TessellatorCommonMethods:

    @PYB11const
    def tessellateQuantized(self,
                            qmesh = "QuantizedTessellation&"):
        """Required for all tessellators:
Compute the quantized tessellation.  This is the basic method all
Tessellator implementations must provide, on which the other tessellation methods
in polytope build."""
        return "void"
