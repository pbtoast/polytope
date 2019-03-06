from mpi4py import *
from polytope import *
from numpy import random

def checkCartesian(mesh, N):
    ncells = mesh.cells.size()
    nnodes = mesh.nodes.size()/2
    nfaces = mesh.faces.size()
    Nc = N*N
    Nn = (N+1)*(N+1)
    Nf = 2*N*(N+1)
    assert ncells == Nc, "Number of cells should be %i, not %i" % (ncells, Nc)
    assert nnodes == Nn, "Number of nodes should be %i, not %i" % (nnodes, Nn)
    assert nfaces == Nf, "Number of faces should be %i, not %i" % (nfaces, Nf)

# ===========================================================================

N = 64
#tess = polytope.TriangleTessellator2d()
tess = polytope.BoostTessellator2d()
degenList = [1.0e-6,
             1.0e-7,
             1.0e-8,
             1.0e-9,
             1.0e-10,
             1.0e-11]

plcPoints = vector_of_double(8)
for i,p in enumerate([0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0]):
    plcPoints[i] = p

plc  = polytope.PLC2d()
plc.facets.resize(4)
for i in range(4):
    plc.facets[i].resize(2)
    plc.facets[i][0] = i
    plc.facets[i][1] = (i+1)%4

points = vector_of_double()
d = 1.0/64
for j in range(N):
    y = (j+0.5)*d
    for i in range(N):
        x = (i+0.5)*d
        points.push_back(x)
        points.push_back(y)


mesh = polytope.Tessellation2d()
for degeneracy in degenList:
    print "Degeneracy = %g" % degeneracy
    tess.tessellateDegenerate(points, plcPoints, plc, degeneracy, mesh)
    checkCartesian(mesh, N)
    mesh.clear()
    print "  PASS"

tess.tessellateDegenerate(points, plcPoints, plc, tess.degeneracy(), mesh)
checkCartesian(mesh, N)
