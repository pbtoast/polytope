from polytope_test_utilities import *
from numpy import *

def positiveID(i):
    if i < 0:  return ~i
    else:      return  i

def orient(p,q):    return sign(p[0]*q[1] - p[1]*q[0])

def orient(x,y,z):  return orient(array(x)-array(y), array(x)-array(z))

def orient(mesh, points, i, j, k):
    assert i < len(list(points))/2
    assert j < mesh.nodes.size()/2
    assert k < mesh.nodes.size()/2
    x = [points[2*i]    , points[2*i+1]]
    y = [mesh.nodes[2*j], mesh.nodes[2*j+1]]
    z = [mesh.nodes[2*k], mesh.nodes[2*k+1]]
    return orient(x,y,z)

def computeSideMesh(points, voro):
    assert points.size()/2 == voro.cells.size()
    sideMesh = polytope.Tessellation2d()
    
    # The nodes are the same. We tack on additional internal nodes (i.e. generators)
    numPointCoords = points.size()
    numNodeCoords  = voro.nodes.size()
    numFaces       = voro.faces.size()
    sideMesh.nodes = voro.nodes
    for p in points:
        sideMesh.nodes.push_back(p)
    assert sideMesh.nodes.size() == numPointCoords + numNodeCoords

    # The number of side mesh cells is #cells * #faces/cell
    numCells = sum([cell.size() for cell in voro.cells])
    sideMesh.cells.resize(numCells)
    
    # The faces are the same. We tack on additional internal faces
    numFaces = voro.faces.size() + numCells
    sideMesh.faces.resize(numFaces)
    sideMesh.faceCells.resize(numFaces)
    for iface,face in enumerate(voro.faces):
        sideMesh.faces[iface] = face

    # Walk the cells and the faces of cells to fill in the new topology
    sideCellIndex = 0
    sideFaceIndex = voro.faces.size()
    for icell, cell in enumerate(voro.cells):
        baseCellIndex = sideCellIndex
        faceOffset = 0
        for face in cell:
            if face < 0:
                iface = ~face
                inode = voro.faces[iface][0]
                sideMesh.faceCells[iface].push_back(~sideCellIndex)
            else:
                iface =  face
                inode = voro.faces[iface][1]
                sideMesh.faceCells[iface].push_back( sideCellIndex)
            nextCellIndex = baseCellIndex + (faceOffset+1)%cell.size()
            assert sideCellIndex < sideMesh.cells.size()
            assert nextCellIndex < sideMesh.cells.size()
            sideMesh.cells[sideCellIndex].push_back(face)
            sideMesh.cells[sideCellIndex].push_back( sideFaceIndex)
            sideMesh.cells[nextCellIndex].push_back(~sideFaceIndex)
            
            assert sideFaceIndex < sideMesh.faces.size(), sideFaceIndex
            assert sideFaceIndex < sideMesh.faceCells.size(), sideFaceIndex
            sideEdge = vector_of_unsigned(2)
            sideEdge[0] = inode
            sideEdge[1] = voro.nodes.size()/2 + icell
            sideMesh.faces[sideFaceIndex] = sideEdge

            sideMesh.faceCells[sideFaceIndex].push_back( sideCellIndex)
            sideMesh.faceCells[sideFaceIndex].push_back(~nextCellIndex)

            sideFaceIndex += 1
            sideCellIndex += 1
            faceOffset += 1

    return sideMesh

# =================================================================================

# Input params
n   = 20
xb  = [0.0, 1.0]
yb  = [0.0, 1.0]
tol = 0.1
smootherType = 1          # 0=small edge nodes, 1=all internal nodes
cycles = 100

# Derived params
Lx = xb[1] - xb[0]
Ly = yb[1] - yb[0]
dx = Lx/sqrt(n)
dy = Ly/sqrt(n)

voroMesh = polytope.Tessellation2d()
T = PyTessellator(serialTessellator = polytope.TriangleTessellator2d(),
                  xlimits = xb,
                  ylimits = yb)

# Construct initial mesh and output it
T.computeRandomTessellation(tessellation = voroMesh, ncells = n)
#T.computeLatticeTessellation(tessellation = voroMesh, nx = 20, ny = 20)
#T.computeHoneycombTessellation(tessellation = voroMesh, nx = 20, ny = 20)
#T.computeRadialTessellation(tessellation = voroMesh, Nrings = 20)
T.outputTessellation(voroMesh, "sideTest", cycle=0)

# Construct the side mesh and output it
sideMesh = computeSideMesh(T.points, voroMesh)
#T.outputTessellation(sideMesh, "sideTest", cycle=1)

polytope.writeTessellation2d(sideMesh, "sideTest", None, None, None, None, cycle=1)
 

