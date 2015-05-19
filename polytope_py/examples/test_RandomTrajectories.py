from polytope_test_utilities import *
from numpy import *
from VelocityFields import *
import time


def newBoundaryPoints(points, dx):
    low = vector_of_double(2)
    low[0] = min([points[2*i  ] for i in range(len(points)/2)])
    low[1] = min([points[2*i+1] for i in range(len(points)/2)])
    #low[0] = min(points)
    #low[1] = low[0]
    plc = polytope.constructConvexHull2d(points, low, dx)
    boundaryPoints = []
    for i,fct in enumerate(plc.facets):
        boundaryPoints.append(points[2*fct[0]  ])
        boundaryPoints.append(points[2*fct[0]+1])
    return boundaryPoints


nx = 20
ny = 20
name = "RandomTrajectories"
dt   = 0.1
Nsteps = 200
robustTessellate = False

mesh = polytope.Tessellation2d()
T = PyTessellator(robustTessellate=robustTessellate)
dx = T.serialTessellator.degeneracy()
N = nx*ny

# Initialize random displacement field
displacements = vectorize_double([0.1*(random.rand()-0.5) for x in range(2*N)])

# Generate an initial lattice tessellation in the unit square
T.computeLatticeTessellation(tessellation=mesh, nx=nx, ny=ny)

# Update the boundary to be the convex hull of the initial points and retessellate
boundaryPoints = newBoundaryPoints(T.points, dx)
T.computeTessellationFromGenerators(mesh, T.points, boundaryPoints=boundaryPoints)
T.outputTessellation(mesh, name, cycle=0)

# Update the background velocity
V = SolidRotation(L = 1.0, v0 = 1.0, origin=[0.5, 0.5])
#V = TaylorGreenVortex(L = 1.0, v0 = 1.0, origin=[0.5, 0.5])

for i in range(Nsteps):
    icycle = i+1
    print "Step %i of %i" % (icycle, Nsteps)
    newPoints = V.update(T.points, dt)
    newPoints = [p + d for (p,d) in zip(newPoints, displacements)]
    boundaryPoints = newBoundaryPoints(vectorize_double(newPoints), dx)
    T.computeTessellationFromGenerators(mesh, newPoints, boundaryPoints=boundaryPoints)
    T.outputTessellation(mesh, name, cycle=icycle)

