from PYB11Generator import *

@PYB11template("int Dimension", "RealType")
class MeshEditor:
    """MeshEditor is used to edit a Tessellation and return a valid topology.
Supported operatations curently include deleting (cells, faces, edges), and 
collapsing degenerate edges by some floating tolerance.
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               mesh = "Tessellation<%(Dimension)s, %(RealType)s>&"):
        "Construct with a tessellation to be edited."
        return

    #...........................................................................
    # Methods
    def deleteCells(self,
                    cellsToDelete = "const std::vector<unsigned>&"):
        "Remote the specified cells (by index), recomputing/ensuring a valid topology in the end."
        return "void"

    def deleteFaces(self,
                    facesToDelete = "const std::vector<unsigned>&"):
        "Remote the specified faces (by index), recomputing/ensuring a valid topology in the end."
        return "void"

    def deleteNodes(self,
                    nodesToDelete = "const std::vector<unsigned>&"):
        "Remote the specified nodes (by index), recomputing/ensuring a valid topology in the end."
        return "void"

    def cleanEdges(self,
                   edgeTol = "const %(RealType)s"):
        "Clean small edges in the mesh based on some tolerance"
        return "void"

    #...........................................................................
    # Attributes
    minEdgesPerFace = PYB11readonly(static=True, doc="Minimum number of edges allowed in a face")

#-------------------------------------------------------------------------------
# Template instantiations
MeshEditor2d = PYB11TemplateClass(MeshEditor, template_parameters=("2", "double"))
MeshEditor3d = PYB11TemplateClass(MeshEditor, template_parameters=("3", "double"))
