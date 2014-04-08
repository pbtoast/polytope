import mpi
from PolytopeModules import *

# Read the generators.
f = open("kuzmin100ascii.txt", "r")
coords = vector_of_double()
print "Reading file..."
for line in f:
    stuff = line.split()
    assert len(stuff) == 3
    assert float(stuff[2]) == 0.0
    for i in xrange(2):
        coords.append(float(stuff[i]))
f.close()
assert len(coords) % 2 == 0
print "Read %i generator positions." % (len(coords) / 2)

# The bounding PLC is defined by the first four generators.
plc = polytope.PLC2d()
plc.facets.resize(4)
for i in xrange(4):
    plc.facets[i].resize(2)
    plc.facets[i][0] = i
    plc.facets[i][1] = (i + 1) % 4
print "Bounding PLC: ", plc

# Build a 2D tessellator.
degen = 0.1
mesh = polytope.Tessellation2d()
tes = polytope.TriangleTessellator2d()
#tes = polytope.BoostTessellator2d()
#tes.tessellate(coords, coords, plc, mesh)
ind = tes.tessellateDegenerate(coords, coords, plc, degen, mesh)
print "Found %i unique generator positions for tessellation." % max(ind)

# Write out a viz file.
polytope.writeTessellation2d(mesh, "Kuzmin_triangle_%g" % degen, None, None, None, None)
