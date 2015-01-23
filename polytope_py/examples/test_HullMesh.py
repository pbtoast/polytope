from mpi4py import *
from polytope_test_utilities import *
from numpy import *

# Params
filename = "generatorData/hullMesh/hullMeshPoints.txt"

# Read the generators.
print "Reading file..."
A = loadtxt(filename)
points = [(x,y) for [x,y] in A]
print "Read %i generator positions." % (len(points))

# Boundary points
boundary = array([0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0])
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
tess = polytope.TriangleTessellator2d()
generators = vector_of_double(2*len(points))
for i,p in enumerate(points):
    generators[2*i  ] = p[0]
    generators[2*i+1] = p[1]

# Tessellate
ind = tess.tessellateDegenerate(generators, plcPoints, plc, tess.degeneracy(), mesh)
print "Found %i unique generator positions." % max(ind)
assert ind.size() == generators.size()/2
#tess.tessellate(generators, plcPoints, plc, mesh)


# Output
cellFields = {}
genx = [generators[2*i  ] for i in range(generators.size()/2)]
geny = [generators[2*i+1] for i in range(generators.size()/2)]
cellFields["gen_x"] = genx
cellFields["gen_y"] = geny
polytope.writeTessellation2d(mesh, "HullMesh", None, None, None, cellFields)
