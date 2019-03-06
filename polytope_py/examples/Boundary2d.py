from numpy import pi, sin, cos, array
from polytope import *
from polytope import vector_of_double as vecd
from polytope import vector_of_vector_of_int as vecveci

# Simply class for handling construction of simple
# boundaries with simple holes. Includes routines for
# testing if a point is located inside a given
# boundary for generator insertion as well as orientation
# checks for outer boundary vertices (CCW) and hole
# vertices (CW)

class Boundary2d:
    def __init__(self):
        self.PLC       = polytope.PLC2d()
        self.PLCpoints = vecd()
        self.center    = None
        self.low       = None
        self.high      = None
        self.outerSet  = False

    # ----------------------------------------------------
    # Worker routine for converting an ordered list of
    # (x,y) positions to STL and PLC data
    # ----------------------------------------------------
    def setOuterBoundary(self, vertices):
        assert len(vertices) % 2 == 0
        assert self.area(vertices) > 0.0
        N = len(vertices)/2
        self.PLCpoints.resize(2*N)
        self.PLC.facets.resize(N)
        xlow = 1.0e80;  xhigh = -1.0e80;
        ylow = 1.0e80;  yhigh = -1.0e80;
        for i in range(N):
            self.PLCpoints[2*i  ] = vertices[2*i  ]
            self.PLCpoints[2*i+1] = vertices[2*i+1]
            xlow  = min(xlow , vertices[2*i  ])
            ylow  = min(ylow , vertices[2*i+1])
            xhigh = max(xhigh, vertices[2*i  ])
            yhigh = max(yhigh, vertices[2*i+1])
            self.PLC.facets[i].push_back(i)
            self.PLC.facets[i].push_back((i+1)%N)
        self.outerSet = True
        self.low  = array([xlow , ylow ])
        self.high = array([xhigh, yhigh])

    # ----------------------------------------------------
    # Worker routine for adding a CW-oriented hole to a
    # geometry, provided the outer boundary and facets
    # are already set.
    # ----------------------------------------------------
    def addHole(self, holeVertices):
        assert self.outerSet
        assert len(holeVertices) % 2 == 0
        assert self.area(holeVertices) < 0.0
        N = len(holeVertices)/2
        assert self.PLCpoints.size() % 2 == 0
        istart = self.PLCpoints.size()/2
        ihole = self.PLC.holes.size()
        self.PLC.holes.push_back(vecveci())
        self.PLC.holes[ihole].resize(N)
        for i in range(N):
            self.PLCpoints.push_back(holeVertices[2*i  ])
            self.PLCpoints.push_back(holeVertices[2*i+1])
            self.PLC.holes[ihole][i].push_back(istart + i      )
            self.PLC.holes[ihole][i].push_back(istart + (i+1)%N)

    # ----------------------------------------------------
    # Some canned hole geometries
    # ----------------------------------------------------

    # Circular hole data
    def addCircularHole(self, radius, center, nfacets): 
        holeVertices = []
        for i in range(nfacets):
            theta = 2.0*pi*(1.0 - i/(nfacets+1.0))
            holeVertices.append(center[0] + radius*cos(theta))
            holeVertices.append(center[1] + radius*sin(theta))
        self.addHole(holeVertices)        

    # Box hole data
    def addBoxHole(self, low, high):
        assert low[0] < high[0] and low[1] < high[1]
        holeVertices = [low [0], low [1], low [0], high[1],
                        high[0], high[1], high[0], low [1]]
        self.addHole(holeVertices)

    # ----------------------------------------------------
    # Some canned geometries
    # ----------------------------------------------------

    # Box outer boundary data
    def initBox(self,
                low = [0.0,0.0],
                high= [1.0,1.0]):
        vertices = [low [0], low [1], high[0], low [1],
                    high[0], high[1], low [0], high[1]]
        self.center = array([0.5*(low[0]+high[0]),
                             0.5*(low[1]+high[1])])
        self.setOuterBoundary(vertices)

    # Circular outer boundary data
    def initCircle(self,
                   radius  = 1.0,
                   center  = [0.0, 0.0],
                   nfacets = 90):
        self.center = array(center)
        vertices = []
        for i in range(nfacets):
            theta = 2.0*pi*i/(nfacets+1.0)
            vertices.append(center[0] + radius*cos(theta))
            vertices.append(center[1] + radius*sin(theta))
        self.setOuterBoundary(vertices)

    # Torus boundary data
    def initTorus(self,
                  outerRadius=1.00,
                  innerRadius=0.25,
                  center=[0.0,0.0],
                  nfacets=90):
        assert outerRadius > innerRadius
        self.initCircle(outerRadius, center, nfacets)
        self.addCircularHole(innerRadius, center, nfacets)

    # M-shaped boundary with two square holes data
    def initMWithHoles(self):
        vertices = [0.0, 0.0,
                    2.0, 0.0,
                    2.0, 2.0,
                    1.0, 1.0,
                    0.0, 2.0]
        self.setOuterBoundary(vertices)
        self.addBoxHole([0.25,0.25], [0.75,0.75])
        self.addBoxHole([1.25,0.25], [1.75,0.75])

    # Unit circle with 5-point star hole boundary data
    def setCircleWithStarHole(self):
        self.initCircle(1.0, [0.0, 0.0], 90)
        theta0 = 2*pi/5
        outerRadius = 0.75
        innerRadius = outerRadius*(sin(theta0/4.0) / sin(3*theta0/4.0))
        holeVertices = []
        for i in range(5):
            theta = pi/2 - i*theta0
            holeVertices.append(outerRadius*cos(theta))
            holeVertices.append(outerRadius*sin(theta))
            theta = pi/2 - i*theta0 - 0.5*theta0
            holeVertices.append(innerRadius*cos(theta))
            holeVertices.append(innerRadius*sin(theta))
        self.addHole(holeVertices)

    # Cardioid boundary data
    def initCardioid(self,
                     z=2,
                     nfacets=90):
        assert z > 0
        vertices = []
        for i in range(nfacets):
            theta = 2.0*pi*i/(nfacets+1)
            vertices.append(z*cos(theta) - cos(2*theta))
            vertices.append(z*sin(theta) - sin(2*theta))
        self.setOuterBoundary(vertices)

    # Trogdor boundary data
    def initTrogdor(self):
        vertices = [5.2, 2.0, 7.0, 1.5, 7.2, 2.7, 8.0, 1.2,
                    8.5, 1.2, 8.5, 3.0, 9.5, 1.6, 10.4, 2.0,
                    9.8, 3.8, 11.4, 3.3, 11.9, 5.6, 11.0, 7.0,
                    9.8, 7.8, 7.0, 9.0, 5.6, 10.5, 4.9, 9.6,
                    3.4, 8.8, 4.0, 8.2, 4.1, 6.5, 5.0, 6.6,
                    5.2, 6.1, 5.0, 6.0, 5.1, 5.5, 5.4, 5.5,
                    5.3, 5.0, 4.0, 5.0, 3.0, 6.0, 1.5, 7.0,
                    1.1, 8.0, 0.7, 8.7, 1.3, 9.6, 0.8, 10.2,
                    0.9, 11.7, 2.0, 12.5, 3.2, 12.3, 3.8, 13.0,
                    4.6, 13.1, 5.2, 14.0, 4.0, 14.3, 3.0, 14.0,
                    2.8, 15.0, 1.9, 15.3, 1.3, 16.2, 0.5, 16.6,
                    2.2, 17.5, 3.6, 15.9, 5.5, 14.3, 6.5, 14.8,
                    7.5, 17.2, 7.9, 19.0, 9.8, 18.5, 9.7, 17.8,
                    9.8, 17.0, 9.2, 16.8, 8.8, 16.0, 8.0, 15.9,
                    7.2, 15.1, 9.6, 15.0, 12.2, 14.0, 14.5, 14.0,
                    14.9, 14.8, 15.5, 14.8, 16.0, 14.2, 15.8, 12.7,
                    12.5, 12.7, 13.1, 11.8, 15.5, 11.0, 15.4, 10.5,
                    13.5, 10.4, 12.0, 11.3, 9.9, 12.1, 8.8, 12.3,
                    8.7, 12.0, 9.0, 11.0, 11.0, 9.2, 14.0, 7.0,
                    14.8, 5.0, 14.0, 3.0, 12.0, 1.1, 10.2, 0.5,
                    8.7, 0.4, 7.2, 0.8]
        self.setOuterBoundary(vertices)

    # ----------------------------------------------------
    # Basic geometric routines
    # ----------------------------------------------------

    # Test if (x,y) lies inside the outer boundary and outside each hole
    def testInside(self, x, y):
        offset = 0
        nsides = self.PLC.facets.size()
        nholes = self.PLC.holes.size()
        isInside = self.inside(x,y,nsides,offset)
        offset += nsides
        for i in range(nholes):
            nsides = self.PLC.holes[i].size()
            isInside ^= self.inside(x,y,nsides,offset)
            offset += nsides
        return isInside

    # Test if an (x,y) point lies inside a PLCpoint collection
    def inside(self, x, y, nsides, offset):
        j = nsides - 1
        isInside = False
        for i in range(nsides):
            ix = 2*(i+offset);  iy = 2*(i+offset)+1
            jx = 2*(j+offset);  jy = 2*(j+offset)+1
            xi = self.PLCpoints[ix];  yi = self.PLCpoints[iy]
            xj = self.PLCpoints[jx];  yj = self.PLCpoints[jy]
            if (((yi <  y and yj >= y) or
                 (yj <  y and yi >= y)) and
                ( xi <= x or  xj <= x)):
                isInside ^= (xi + (y-yi)/(yj-yi)*(xj-xi) < x)
            j = i
        return isInside

    # Compute the area of a 2D polygon
    def area(self, vertices):
        assert len(vertices) % 2 == 0
        N = len(vertices)/2
        area = 0.0
        for ii in range(N):
            i0 = ii;  i1 = (i0+1) % N
            x0 = vertices[2*i0];  y0 = vertices[2*i0+1]
            x1 = vertices[2*i1];  y1 = vertices[2*i1+1]
            area += (x0*y1 - y0*x1)
        return area
