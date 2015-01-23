from polytope_test_utilities import *
from numpy import *
import time

# Params
N = 9
xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 1.0

# Derived
xlimits = array([xmin, xmax])
ylimits = array([ymin, ymax])

# Initialize Tessellator stuff
partitioner = VoronoiPartitioner(xlimits = xlimits,
                                 ylimits = ylimits,
                                 lloydIterations = 10)
mesh = polytope.Tessellation2d()
T = PyTessellator(serialTessellator = polytope.BoostTessellator2d(),
                  xlimits           = xlimits,
                  ylimits           = ylimits,
                  partitioner       = partitioner)

# Tessellate!
for i in range(N):
    print "%i of %i" % ((i+1), N)
    nx = 2**(i+1)
    ny = 2**(i+1)
    mesh.clear()
    T.computeHoneycombTessellation(mesh, nx, ny)
    T.outputTessellation(mesh, "Honeycomb", cycle=i)
