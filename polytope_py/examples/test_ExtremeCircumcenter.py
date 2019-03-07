from polytope_test_utilities import *
from mpi4py import *
from polytope import *

def outputTessellation(tessellation, points, name,
                       nodeFields = None,
                       edgeFields = None,
                       faceFields = None,
                       cellFields = None,
                       cycle = 0):
    if cellFields is None:
        cellFields = {}
    genx = [points[2*i  ] for i in range(points.size()/2)]
    geny = [points[2*i+1] for i in range(points.size()/2)]
    cellFields["gen_x"] = genx
    cellFields["gen_y"] = geny
    polytope.writeTessellation2d(tessellation, name, nodeFields, edgeFields, faceFields, cellFields, 
                                 cycle=cycle)

def circ(A,B,C):
    D = 2.0*(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1])) + 1.0e-80
    X = ((A[0]*A[0] + A[1]*A[1])*(B[1]-C[1]) + 
         (B[0]*B[0] + B[1]*B[1])*(C[1]-A[1]) + 
         (C[0]*C[0] + C[1]*C[1])*(A[1]-B[1]))/D;
    Y = ((A[0]*A[0] + A[1]*A[1])*(C[0]-B[0]) + 
         (B[0]*B[0] + B[1]*B[1])*(A[0]-C[0]) + 
         (C[0]*C[0] + C[1]*C[1])*(B[0]-A[0]))/D;
    return X,Y

L = 1.0
A = 8
staticPts = [0.0,  0.5*L,
             0.0, -0.5*L,
              -L,  0.5*L,
              -L, -0.5*L,
             0.0,  0.0001*L,
             0.0, -0.0001*L]
xmax = (8*A - 0.5)*L
ymax =  8*A*L

tessellator = polytope.TriangleTessellator2d()
X = -L
for i in range(40):
    X += (L/2)/(2**i)
    points = vector_of_double()
    for p in staticPts:  points.push_back(p)
    points.push_back(X)
    points.push_back(0.0)
    mesh = polytope.Tessellation2d()
    tessellator.tessellate(points, mesh)
    outputTessellation(mesh, points, "ExtremeCircumcenter", cycle=i)
    (xc1,yc1) = circ([X,0], [0.0, 0.5*L], [0.0, -0.5*L])
    (xc2,yc2) = circ([X,0], [xmax, ymax], [0.0, 0.5*L])
    print "%2i %8.8f %8.8f" % (i, xc1, xc2)
    #if xc1 < xc2:  print "%2i %8.8f  " % (i, xc1)
    #else:          print "%2i %8.8f !" % (i, xc2)
    
