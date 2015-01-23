from polytope_test_utilities import *
from numpy import *

# Params
radius             = 0.5
center             = array([0.5, 0.5])
Nrings             = 160
centralGenerator   = True
nfold              = 6
boundaryGenerators = False
printResults       = False

# Derived
xlimits = center + radius*array([-1,1])
ylimits = center + radius*array([-1,1])
nArcs   = nfold*(Nrings-1)
boundaryPoints = []
for i in range(nArcs):
    theta = 2*pi*i/nArcs
    boundaryPoints.append(center[0] + radius*cos(theta))
    boundaryPoints.append(center[1] + radius*sin(theta))

# Initialize Tessellator stuff
partitioner = VoronoiPartitioner(boundaryPoints = boundaryPoints,
                                 xlimits = xlimits,
                                 ylimits = ylimits,
                                 lloydIterations = 20)
mesh = polytope.Tessellation2d()
T = PyTessellator(serialTessellator = polytope.BoostTessellator2d(),
                  boundaryPoints    = boundaryPoints,
                  xlimits           = xlimits,
                  ylimits           = ylimits,
                  partitioner       = partitioner,
                  timer             = True,
                  boundaryGenerators= boundaryGenerators)

# Tessellate!
T.computeRadialTessellation(tessellation       = mesh,
                            Nrings             = Nrings,
                            center             = center,
                            centralGenerator   = centralGenerator,
                            nfold              = nfold,
                            smallEdgeTolerance = 0.0)

polytope.writeTessellation2d(mesh, "RadialTessellation", None, None, None, None)

if printResults:
   x = [T.points[2*i  ] for i in range(T.points.size()/2)]
   y = [T.points[2*i+1] for i in range(T.points.size()/2)]
   fileName = "radialGenerators_%irings_%iprocs-%i" % (Nrings, MPI.COMM_WORLD.size, MPI.COMM_WORLD.rank)
   f = open(fileName, "w")
   f.write(("#" + 2*"%20s   " + "\n") % ("x", "y"))
   for xi,yi in zip(x,y):
      f.write((" " + 2*"%20.16f   " + "\n") % (xi,yi))
   f.close()

