from mpi4py import *
from polytope_test_utilities import *
from numpy import *
comm = MPI.COMM_WORLD

assert comm.Get_size() == 20

# Params
numGenerators = 100000
filename = "generatorData/impact/impactPoints_%i.txt" % comm.Get_rank()

# Read the generators.
print "Reading file..."
A = loadtxt(filename)
points = [(x,y) for [x,y] in A]
print "Read %i generator positions." % (len(points))

# Boundary points
boundary = array([0.0, 0.0, 2.0, 0.0, 2.0, 2.0, 0.0, 2.0])
plc = polytope.PLC2d()
plc.facets.resize(4)
for i in xrange(4):
    plc.facets[i].resize(2)
    plc.facets[i][0] = i
    plc.facets[i][1] = (i + 1) % 4
plcPoints = vector_of_double(8)
for i,p in enumerate(boundary):
    plcPoints[i] = p


# Initialize the tessellator
mesh = polytope.Tessellation2d()
triangle = polytope.TriangleTessellator2d()
tess = polytope.DistributedTessellator2d(triangle, False, True)
generators = vector_of_double(2*len(points))
for i,p in enumerate(points):
    generators[2*i  ] = p[0]
    generators[2*i+1] = p[1]

# Tessellate
tol = 1.0e-6
tess.tessellate(generators, plcPoints, plc, mesh)
polytope.writeTessellation2d(mesh, "Impact", None, None, None, None)
