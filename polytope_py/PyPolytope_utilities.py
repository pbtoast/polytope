from numpy import *
from mpi4py import *
from polytope import *

#-------------------------------------------------------------------------------
# Convert a list to Pybindgen-wrapped std vector
#-------------------------------------------------------------------------------
def vectorize(inputList):
    if all(isinstance(x,int) for x in inputList):
        return vectorize_int(inputList)
    else:
        return vectorize_double(inputList)

#-------------------------------------------------------------------------------
# Vectorize for integer types
#-------------------------------------------------------------------------------
def vectorize_int(inputList):
    result = vector_of_int()
    result.resize(len(inputList))
    for i,p in enumerate(inputList):
        result[i] = p
    return result

#-------------------------------------------------------------------------------
# Vectorize for float types
#-------------------------------------------------------------------------------
def vectorize_double(inputList):
    result = vector_of_double()
    result.resize(len(inputList))
    for i,p in enumerate(inputList):
        result[i] = p
    return result

#-------------------------------------------------------------------------------
# Turn set of [xmin,xmax], [ymin,ymax] into CCW ordered boundary points
#-------------------------------------------------------------------------------
def boxListFromLimits(xlimits, ylimits):
    assert len(xlimits) == len(ylimits)
    return boxList(xlimits[0], xlimits[1], ylimits[0], ylimits[1])

#-------------------------------------------------------------------------------
# Turn set of min/max values into CCW ordered (x1,y1,x2,y2,...) box boundary
#-------------------------------------------------------------------------------
def boxList(xmin, xmax, ymin, ymax):
    assert (xmax > xmin and ymax > ymin)
    return [xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax]

#-------------------------------------------------------------------------------
# Compute points on a lattice
#-------------------------------------------------------------------------------
def computeLatticePoints(xmin, ymin, nx, ny, dx, dy):
    xps = [xmin + (i+0.5)*dx for i in range(nx)]
    yps = [ymin + (i+0.5)*dy for i in range(ny)]
    pps = []
    for y in yps:
        for x in xps:
            pps.append([x,y])
    return pps

#-------------------------------------------------------------------------------
# Compute squared distance of two vectors
#-------------------------------------------------------------------------------
def distance2(a,b):
    assert len(a) > 0
    assert len(a) == len(b)
    d = 0.0
    for ai,bi in zip(a,b):
        d += (ai-bi)**2
    return d

#-------------------------------------------------------------------------------
# Return a PLC from an ordered list of points [x1,y1,x2,y2,...]
#-------------------------------------------------------------------------------
def computePLCfromList2d(PLCpoints):
    assert len(PLCpoints) % 2 == 0
    N = len(PLCpoints)/2
    PLC = polytope.PLC2d()
    PLC.facets.resize(N)
    for i in range(N):
        PLC.facets[i].resize(2)
        PLC.facets[i][0] = i
        PLC.facets[i][1] = (i+1) % N
    return PLC

#-----------------------------------------------------------------------------
# Perform nStep iterations of Lloyd's algorithm to smooth bounded tessellation
#-----------------------------------------------------------------------------
def smoothTessellation(tessellation, 
                       tessellator,
                       points,
                       PLCpoints,
                       PLC,
                       iterations  = 10):
    assert not PLC.empty()
    assert points.size() > 0
    assert PLCpoints.size() > 0
    for i in range(iterations):
        tessellation.clear()
        tessellator.tessellate(points, PLCpoints, PLC, tessellation)
        centroids = computeCentroids(tessellation)
        assert len(centroids) == points.size()
        for j,c in enumerate(centroids):
            points[j] = 0.5*(points[j] + c)
        #polytope.writeTessellation2d(tessellation, "DEBUG_smoothTessellation",
        #                             None, None, None, None, cycle=i)
    tessellation.clear()
    tessellator.tessellate(points, PLCpoints, PLC, tessellation)
    
#-----------------------------------------------------------------------------
# Compute centroids of Polytope tessellation cells
#-----------------------------------------------------------------------------
def computeCentroids(tessellation):
    centroids = []
    for icell,cell in enumerate(tessellation.cells):
        area = 0.0
        cx = 0.0
        cy = 0.0
        for ftmp in cell:
            if ftmp < 0:  
                iface = ~ftmp
                n0 = tessellation.faces[iface][1]
                n1 = tessellation.faces[iface][0]
            else:         
                iface = ftmp
                n0 = tessellation.faces[iface][0]
                n1 = tessellation.faces[iface][1]
            x0 = tessellation.nodes[2*n0  ]
            y0 = tessellation.nodes[2*n0+1]
            x1 = tessellation.nodes[2*n1  ]
            y1 = tessellation.nodes[2*n1+1]
            d = x0*y1 - y0*x1
            area += d
            cx   += d*(x0+x1)
            cy   += d*(y0+y1)
            #assert area > 0.0, "Negative area computed in centroid calculation"
        area /= 2.0
        cx   /= (6*area)
        cy   /= (6*area)
        centroids.append(cx)
        centroids.append(cy)
    return centroids
