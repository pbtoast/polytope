from numpy import *
from polytope_test_utilities import *
import time

# ==============================================================
# Weight Functions
# Weights are computed using exp(-100*(F(x,y)))
# Generators are dragged towards the F(x,y)=0 level curve
# ==============================================================
def Constant(x,y):  return 1.0
def Point(x,y):     return (x-0.5)**2 + (y-0.5)**2
def Ring(x,y):      return abs((x-0.5)**2 + (y-0.5)**2 - 0.125)
def Diagonal(x,y):  return (y-x)**2
def Square(x,y):    return abs(max(abs(x-0.5), abs(y-0.5)) - 0.25)
def Hyperbola(x,y): return abs((y-0.5)**2 - (x-0.5)**2)

# ---------------------------------------------------------------

def computeWeight(x,y,Function):
    return exp(-100*( Function(x,y) ))

def computeBarycenters(tessellation):
    xb       = [0.0 for i in range(tessellation.cells.size())]
    yb       = [0.0 for i in range(tessellation.cells.size())]
    numfaces = [cell.size() for cell in tessellation.cells]
    for iface,face in enumerate(tessellation.faces):
       fc = tessellation.faceCells[iface]
       xf = 0.5*(tessellation.nodes[2*face[0]  ] + tessellation.nodes[2*face[1]  ])
       yf = 0.5*(tessellation.nodes[2*face[0]+1] + tessellation.nodes[2*face[1]+1])
       for c in fc:
          if c < 0: icell = ~c
          else:     icell = c
          xb[icell] += xf
          yb[icell] += yf
    return [(x/n, y/n) for x,y,n in zip(xb,yb,numfaces)]
       
def computeWeightedCentroids(tessellation, weightFunction):
    barycenters = computeBarycenters(tessellation)
    centroids = []
    for icell, cell in enumerate(tessellation.cells):
        area = 0.0
        xc = 0.0
        yc = 0.0
        for ftmp in cell:
            if ftmp < 0:  
                n0 = tessellation.faces[~ftmp][1]
                n1 = tessellation.faces[~ftmp][0]
            else:         
                n0 = tessellation.faces[ ftmp][0]
                n1 = tessellation.faces[ ftmp][1]
            x0 = tessellation.nodes[2*n0  ]
            y0 = tessellation.nodes[2*n0+1]
            x1 = tessellation.nodes[2*n1  ]
            y1 = tessellation.nodes[2*n1+1]
            xb = barycenters[icell][0]
            yb = barycenters[icell][1]
            xt = (x0 + x1 + xb)/3.0
            yt = (y0 + y1 + yb)/3.0
            xe = (x0 + x1)/2.0
            ye = (y0 + y1)/2.0
            w = computeWeight(xe, ye, weightFunction)
            d  = ((x0-xb)*(y1-yb) - (x1-xb)*(y0-yb))*w
            area += d
            xc   += d*xt
            yc   += d*yt
        area /= 2.0
        xc   /= (2*area)
        yc   /= (2*area)
        centroids.append(xc)
        centroids.append(yc)
    return centroids

# ------------------------------------------------------------

# Params
numIterations = 10
nx            = 100
ny            = 100
Function      = Square     #Constant, Point, Ring, Diagonal, Square, Hyperbola, 
w             = 1.0        #Weight factor on centroid over orig. location
dumpCycle     = 1
serialTess    = "Triangle" #Triangle, Boost

# Derived params
numGenerators = nx*ny

# Checks
assert 0.0 <= w <= 1.0
assert Function in (Constant, Point, Ring, Diagonal, Square, Hyperbola)
assert serialTess in ("Triangle", "Boost")

# Initialize the partitioner
mesh = polytope.Tessellation2d()
T = PyTessellator(serialTessellator = serialTess,
                  xlimits           = [0.0, 1.0],
                  ylimits           = [0.0, 1.0])

# Initial mesh
T.computeLatticeTessellation(tessellation = mesh, nx=nx, ny=ny)
T.outputTessellation(mesh, "WeighedCentroid", cycle=0)

# Cycle:
for k in range(numIterations):
   icycle = k+1
   if MPI.COMM_WORLD.rank == 0: print "%i of %i" % (icycle, numIterations)
   centroids = computeWeightedCentroids(mesh, Function)
   newPoints = [w*c + (1-w)*p for i,(p,c) in enumerate(zip(list(T.points), centroids))]
   T.computeTessellationFromGenerators(mesh, newPoints)
   if (icycle % dumpCycle == 0):  T.outputTessellation(mesh, "WeighedCentroid", cycle=icycle)

