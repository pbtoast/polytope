from polytope_test_utilities import *
from numpy import *
import time

def outputTessellation(tessellation, name, points, cycle = 0):
    cellFields = {}
    genx = [points[2*i  ] for i in range(len(list(points))/2)]
    geny = [points[2*i+1] for i in range(len(list(points))/2)]
    cellFields["gen_x"] = genx
    cellFields["gen_y"] = geny
    polytope.writeTessellation2d(tessellation, name, None, None, None, cellFields, cycle = cycle)



N = 20
Ntests = 20
xmin = 0.0
ymin = 0.0
name = "ConvexHull"

tessellator = polytope.TriangleTessellator2d()
dx = tessellator.degeneracy()

for itest in range(Ntests):
    points = vector_of_double(2*N)
    for i in range(2*N):
        points[i] = random.rand()

    lowerBounds = vector_of_double(2)
    lowerBounds[0] = xmin
    lowerBounds[1] = ymin
    
    plc = polytope.constructConvexHull2d(points, lowerBounds, dx)

    nVerts = plc.facets.size()
    plcPoints = vector_of_double(2*nVerts)
    for i,fct in enumerate(plc.facets):
        plcPoints[2*i  ] = points[2*fct[0]  ]
        plcPoints[2*i+1] = points[2*fct[0]+1]
        fct[0] = i
        fct[1] = (i+1) % nVerts
        
    mesh = polytope.Tessellation2d()
    tessellator.tessellate(points, plcPoints, plc, mesh)
    outputTessellation(mesh, name, points, itest)
