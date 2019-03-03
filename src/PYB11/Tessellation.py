from PYB11Generator import *

@PYB11template("Dimension", "RealType")
class Tessellation:
    "A basic descriptor class for a topologically-consistent arbitrary poly(gonal/hedral) mesh."

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Methods
    def clear(self):
        "Clears the tessellation, emptying it of all data."
        return "void"

    @PYB11const
    def empty(self):
        """Returns true if the tessellation is empty (not defined), 
false otherwise."""
        return "bool"

    def computeNodeCells(self):
        "Find the set of cells that touch each mesh node."
        return "std::vector<std::set<unsigned>>"

    def computeCellToNodes(self):
        "Collect the nodes around each cell"
        return "std::vector<std::set<unsigned>>"

    @PYB11implementation("[](const Tessellation<%(Dimension)s, %(RealType)s>& self) { std::stringstream ss; ss << self; return ss.str();}")
    def __str__(self):
        return "std::string"

    #...........................................................................
    # Attributes
    nodes = PYB11readwrite(doc="""An array of (Dimension*numNodes) values containing components of 
node positions. The components are stored in node-major order and 
the 0th component of the ith node appears in nodes[Dimension*i].""")

    cells = PYB11readwrite(doc="""This two-dimensional array defines the cell-face topology of the 
mesh. A cell has an arbitrary number of faces in 2D and 3D.
cells[i][j] gives the index of the jth face of the ith cell.
A negative face index indicates the actual face index is the 1's 
complement of the value (~cells[i][j]) and the face is oriented
with an inward pointing normal for cells[i].""")

    faces = PYB11readwrite(doc="""This two-dimensional array defines the topology of the faces of the 
mesh. A face has an arbitrary number of nodes in 3D and 2 nodes in 2D. 
faces[i][j] gives the index of the jth node of the ith face.
Nodes for a given face are arranged counterclockwise around the face
viewed from the "positive" (outside) direction.""")

    infNodes = PYB11readwrite(doc="""Indices of all nodes in an unbounded tessellation that are
"infinite." An infinite node is the termination point on a spherical
surface of a ray (tessellation edge) going out to infinity.""")

    infFaces = PYB11readwrite(doc="""Array of face indices: 0 for interior faces and 1 for faces at
"infinity" for unbounded tessellations. The infinite face connects
the collection of infinite nodes for a given unbounded cell.""")

    faceCells = PYB11readwrite(doc="""An array of cell indices for each face, i.e., the cells that share
the face.
For a given cell there will be either 1 or 2 cells -- the cases with 1
cell indicate a face on a boundary of the tessellation.""")

    convexHull = PYB11readwrite(doc="""A PLC connecting the generating points belonging to the convex hull 
of the point distribution. Not all Tessellators hand back the convex 
hull, so this may be empty, in which case you must compute the convex 
hull yourself.""")

    neighborDomains = PYB11readwrite(doc="""Parallel data structure: the set of neighbor domains this portion of
the tessellation is in contact with.""")

    sharedNodes = PYB11readwrite(doc="""Parallel data structure: the nodes this domain shares with
each neighbor domain.
NOTE: we implicitly assume that any domains of rank less than ours we
      are receiving from, while any domains of greater rank we send
      to.""")

    sharedFaces = PYB11readwrite(doc="""Parallel data structure: the faces this domain shares with
each neighbor domain.
NOTE: we implicitly assume that any domains of rank less than ours we
      are receiving from, while any domains of greater rank we send
      to.""")

#-------------------------------------------------------------------------------
# Template instantiations
Tessellation2d = PYB11TemplateClass(Tessellation, template_parameters=("2", "double"))
Tessellation3d = PYB11TemplateClass(Tessellation, template_parameters=("3", "double"))
