from numpy import *
from polytope_test_utilities import *
from VelocityFields import *

def positiveIndex(i):
    if i < 0:  return ~i
    else:      return  i

def computeTGVvelocities(generators):
    velocities = []
    for i in range(len(list(generators))/2):
        u =  0.5*(cos(2*pi*(generators[2*i]-0.25))*sin(2*pi*(generators[2*i+1]-0.25)))
        v = -0.5*(sin(2*pi*(generators[2*i]-0.25))*cos(2*pi*(generators[2*i+1]-0.25)))
        velocities.append(u)
        velocities.append(v)
    return velocities

def moveGenerators(generators, dt):
    velocities = computeTGVvelocities(generators)
    halfTimeGenerators = [g + 0.5*dt*v for (g,v) in zip(generators, velocities)]
    velocities = computeTGVvelocities(halfTimeGenerators)
    return [g + dt*v for (g,v) in zip(generators, velocities)]

def computeBarycenters(tessellation):
    barycenters = []
    for cell in tessellation.cells:
        cx = 0.0
        cy = 0.0
        N  = 0
        for ftmp in cell:
            if ftmp < 0:  n0 = tessellation.faces[~ftmp][1]
            else:         n0 = tessellation.faces[ ftmp][0]
            x0 = tessellation.nodes[2*n0  ]
            y0 = tessellation.nodes[2*n0+1]
            cx += x0
            cy += y0
            N  += 1
        assert N > 0
        cx /= N
        cy /= N
        barycenters.append(cx)
        barycenters.append(cy)
    return barycenters


def computeVoronoiDual(points, voro):
    assert points.size()/2 == voro.cells.size()
    dual = polytope.Tessellation2d()
    
    # The nodes are the generating points
    dual.nodes.resize(points.size())
    for i,p in enumerate(points):
        dual.nodes[i] = p

    # Pull out the boundary nodes. They won't have dual cells
    boundaryNodes = []
    for iface,face in enumerate(voro.faces):
        assert face.size() == 2
        inode1 = face[0]
        inode2 = face[1]
        if voro.faceCells[iface].size() == 1:
            boundaryNodes.append(inode1)
            boundaryNodes.append(inode2)
    boundaryNodes = unique(boundaryNodes)
    numTriangles = voro.nodes.size()/2 - len(boundaryNodes)
    assert numTriangles > 0
    dual.cells.resize(numTriangles)

    # Map each internal voronoi node to a dual cell index
    voroNodeToDualCell = {}
    icell = 0
    for i in range(voro.nodes.size()/2):
        if i not in boundaryNodes:
            voroNodeToDualCell[i] = icell
            icell += 1
    
    # Walk the faces and fill in the rest of the topology
    dualFaceIndex = 0
    for iface,face in enumerate(voro.faces):
        assert face.size() == 2
        inode1 = face[0]
        inode2 = face[1]
        faceCells = voro.faceCells[iface]
        if faceCells.size() == 2:
            if faceCells[0] < 0:
                assert faceCells[1] >= 0
                icell1 = ~faceCells[0]
                icell2 =  faceCells[1]
            else:
                assert faceCells[0] >= 0
                icell1 = ~faceCells[1]
                icell2 =  faceCells[0]
            dualEdge = vector_of_unsigned(2)
            dualEdge[0] = icell1
            dualEdge[1] = icell2
            dual.faces.push_back(dualEdge)
            dualEdgeCells = vector_of_int()
            if inode1 not in boundaryNodes:
                index = voroNodeToDualCell[inode1]
                dual.cells[index].push_back( dualFaceIndex)
                dualEdgeCells.push_back(index)
            if inode2 not in boundaryNodes:
                index = voroNodeToDualCell[inode2]
                dual.cells[index].push_back(~dualFaceIndex)
                dualEdgeCells.push_back(~index)
            assert dualEdgeCells.size() > 0
            dual.faceCells.push_back(dualEdgeCells)
            dualFaceIndex += 1
    return dual, voroNodeToDualCell

def dualSmoother(generators, mesh, smootherType=0, tol=0.1):
    # Walk the Voronoi mesh and compute min/max/avg edge lengths
    minLength =  1.0e99
    maxLength = -1.0e99
    avgLength =  0.0
    edgeLengths = []
    boundaryNodes = []
    for iface,face in enumerate(mesh.faces):
        inode1 = face[0]
        inode2 = face[1]
        n1 = [mesh.nodes[2*inode1], mesh.nodes[2*inode1+1]]
        n2 = [mesh.nodes[2*inode2], mesh.nodes[2*inode2+1]]
        edgeLength = sqrt((n1[0]-n2[0])**2 + (n1[1]-n2[1])**2)
        minLength = min(minLength, edgeLength)
        maxLength = max(maxLength, edgeLength)
        avgLength += edgeLength
        edgeLengths.append(edgeLength)
        if mesh.faceCells[iface].size() == 1:
            boundaryNodes.append(inode1)
            boundaryNodes.append(inode2)
    avgLength /= mesh.faces.size()
    boundaryNodes = unique(boundaryNodes)
    internalNodes = setdiff1d(arange(mesh.nodes.size()/2), boundaryNodes)

    # Flag all the small edges
    flaggedNodes = []
    for i,L in enumerate(edgeLengths):
        if L < tol*avgLength:
            flaggedNodes.append(mesh.faces[i][0])
            flaggedNodes.append(mesh.faces[i][1])
    flaggedNodes = unique(flaggedNodes)

    # Construct dual mesh and output it
    dualMesh, voroNodeToDualCell = computeVoronoiDual(generators, mesh)
    #polytope.writeTessellation2d(dualMesh, "dualMesh", None, None, None, None)

    # Compute the barycenters of the dual mesh cells
    dualBarycenters = computeBarycenters(dualMesh)

    if smootherType == 0:
        nodeList = flaggedNodes
    elif smootherType == 1:
        nodeList = internalNodes
    else:
        assert False

    # Now move the flagged nodes to "more optimal" positions
    for inode in nodeList:
        if inode in voroNodeToDualCell:
            dualCellIndex = voroNodeToDualCell[inode]
            assert dualCellIndex < dualMesh.cells.size()
            assert dualCellIndex < len(dualBarycenters)/2
            mesh.nodes[2*inode  ] = dualBarycenters[2*dualCellIndex  ] 
            mesh.nodes[2*inode+1] = dualBarycenters[2*dualCellIndex+1]
            
    return

# =================================================================================

# Input params
n   = 100
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
#T.computeRandomTessellation(tessellation = voroMesh, ncells = n)
T.computeLatticeTessellation(tessellation = voroMesh, nx = 20, ny = 20)
#T.computeHoneycombTessellation(tessellation = voroMesh, nx = 20, ny = 20)
#T.computeRadialTessellation(tessellation = voroMesh, Nrings = 20)
T.outputTessellation(voroMesh, "voroMesh_unsmoothed", cycle=0)
T.outputTessellation(voroMesh, "voroMesh_smoothed",   cycle=0)
 
dt = sqrt(2)*min(dx,dy)/2.0
for icycle in range(cycles):
    newPoints = moveGenerators(T.points, dt)
    T.computeTessellationFromGenerators(tessellation = voroMesh, points = newPoints)
    T.outputTessellation(voroMesh, "voroMesh_unsmoothed", cycle=icycle+1)
    dualSmoother(T.points, voroMesh, smootherType=smootherType, tol=tol)
    T.outputTessellation(voroMesh, "voroMesh_smoothed",   cycle=icycle+1)


