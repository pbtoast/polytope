from polytope_test_utilities import *
from numpy import *

# Params
numGenerators = 100000
filename = "generatorData/kuzmin/kuzmin100ascii.txt"

# Read the generators.
print "Reading file..."
A = loadtxt(filename)
points = [(x,y) for [x,y,z] in A[0:numGenerators,:]]
assert len(points) == numGenerators
print "Read %i generator positions." % (len(points))

# Boundary points
boundary = 5000.0*array([-1, -1, 1, -1, 1, 1, -1, 1])
plc = polytope.PLC2d()
plc.facets.resize(4)
for i in xrange(4):
    plc.facets[i].resize(2)
    plc.facets[i][0] = i
    plc.facets[i][1] = (i + 1) % 4

# Initialize the tessellator
mesh = polytope.Tessellation2d()
tess = polytope.TriangleTessellator2d()
generators = vector_of_double(2*len(points))
for i,p in enumerate(points):
    generators[2*i  ] = p[0]
    generators[2*i+1] = p[1]

# Tessellate
tol = 1.0e-6
ind = tess.tessellateDegenerate(generators, tol, mesh)
#ind = tess.tessellateDegenerate(generators, generators, plc, tol, mesh)
print "Found %i unique generator positions." % max(ind)
polytope.writeTessellation2d(mesh, "Kuzmin", None, None, None, None)
