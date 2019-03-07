from numpy import *
from mpi4py import *
from PyPolytope_utilities import *
from PyPartitioner import *
from polytope import *
import time

#================================================================================
# PyTessellator
# Class to encapsulate partitioning points, handling boundaries, and tessellating
#================================================================================
class PyTessellator:
    def __init__(self,
                 serialTessellator = None,
                 useSerial         = False,
                 boundaryPoints    = None,
                 xlimits           = [0.0, 1.0],
                 ylimits           = [0.0, 1.0],
                 partitioner       = None,
                 timer             = False,
                 boundaryGenerators= False,
                 unbounded         = False,
                 robustTessellate  = True):

        # Store the serial tessellator - default is Triangle
        if serialTessellator is not None:
            if   serialTessellator == "Boost":
                self.serialTessellator = polytope.BoostTessellator2d()
            elif serialTessellator == "Triangle":
                self.serialTessellator = polytope.TriangleTessellator2d()
            else:
                self.serialTessellator = serialTessellator
        else:
            self.serialTessellator = polytope.TriangleTessellator2d()

        # Construct the distributed tessellator that actually gets used.
        # Default is Distributed
        if useSerial:
            self.tessellator = polytope.SerialDistributedTessellator2d(self.serialTessellator, False, True)
        else:
            self.tessellator = polytope.DistributedTessellator2d(self.serialTessellator, False, True)

        # Store the domain limits
        self.xlimits = xlimits
        self.ylimits = ylimits
        
        # Store the input boundary as vector_of_double PLCpoints and a faceted PLC
        # Default is to initialize x- and y-limits
        if boundaryPoints is None:
            boundaryPoints = boxListFromLimits(xlimits, ylimits)
        self.boundaryPoints = boundaryPoints
        self.PLC            = computePLCfromList2d(boundaryPoints)
        self.PLCpoints      = vectorize_double(boundaryPoints)
        
        # Initialize generators
        self.points = vector_of_double()

        # MPI bidness
        self.comm = MPI.COMM_WORLD

        # Time the tesselate calls and return results
        self.timer = timer

        # Are the boundary points also generators?
        self.boundaryGenerators = boundaryGenerators

        # Use the degeneracy-safe (robust) tessellate call?
        self.robustTessellate = robustTessellate

        # Partitioner - default is Voronoi partitioner
        if partitioner is None:
            self.partitioner = VoronoiPartitioner(self.serialTessellator,
                                                  boundaryPoints,
                                                  xlimits,
                                                  ylimits)
        else:
            self.partitioner = partitioner
            
        # Unbounded
        self.unbounded = unbounded

        return

    #-------------------------------------------------------------------------------
    # Set the degeneracy threshold for the serial tessellator
    #-------------------------------------------------------------------------------
    def degeneracy(self, val):
        self.serialTessellator.degeneracy(val)

    #-------------------------------------------------------------------------------
    # Generic tessellate routine
    #-------------------------------------------------------------------------------
    def tessellate(self, tessellation):
        if self.boundaryGenerators:
            allPoints = vectorize_double(list(self.points)+list(self.PLCpoints))
        else:
            allPoints = self.points

        if not tessellation.empty():
            tessellation.clear()

        if self.timer:  t0 = time.clock()
        if self.unbounded:
            self.tessellator.tessellate(self.points, tessellation)
        else:
            if self.robustTessellate:
                self.tessellator.tessellateDegenerate(allPoints, 
                                                      self.PLCpoints, 
                                                      self.PLC, 
                                                      self.tessellator.degeneracy(),
                                                      tessellation)
            else:
                #self.tessellator.tessellate(allPoints, self.PLCpoints, self.PLC, tessellation)
                self.tessellator.tessellateNormalized(allPoints, self.PLCpoints, self.PLC, tessellation)



        if self.timer:
            t1 = time.clock()
            nPoints = self.comm.allreduce(self.points.size()/2, op=MPI.SUM)
            if self.comm.rank == 0:
                print "%g generators required %g seconds" % \
                    (nPoints,t1-t0)
    
    #-------------------------------------------------------------------------------
    # Lattice tessellation utility function
    # Generate points on a lattice, partition and store them, and tessellate
    #-------------------------------------------------------------------------------
    def computeLatticeTessellation(self, tessellation, nx, ny):
        # Lattice generators
        dx = (self.xlimits[1] - self.xlimits[0])/nx
        dy = (self.ylimits[1] - self.ylimits[0])/ny
        pps = computeLatticePoints(self.xlimits[0], self.ylimits[0], nx, ny, dx, dy)

        # Partition the points
        self.points = self.partitioner.computePartition(pps)

        # Make sure the boundary data corresponds to the specified limits
        boundaryPoints = boxListFromLimits(self.xlimits, self.ylimits)
        self.boundaryPoints = boundaryPoints
        self.PLC            = computePLCfromList2d(boundaryPoints)
        self.PLCpoints      = vectorize_double(boundaryPoints)

        # Tessellate!
        self.tessellate(tessellation)
                                                                      

    #-------------------------------------------------------------------------------
    # Random tessellation utility function
    # Generate random points in a box, partition and store them, and tessellate
    #-------------------------------------------------------------------------------
    def computeRandomTessellation(self, tessellation, ncells, 
                                  seed = 13927391):
        # Random generators
        random.seed(seed)
        pps = []
        for i in range(ncells):
            xrand = self.xlimits[0] + (self.xlimits[1]-self.xlimits[0])*random.rand()
            yrand = self.ylimits[0] + (self.ylimits[1]-self.ylimits[0])*random.rand()
            pps.append([xrand,yrand])
            
        # Partition the generators
        self.points = self.partitioner.computePartition(pps)
                
        # Make sure the boundary data is non-empty
        assert not self.PLC.empty()
        assert self.PLCpoints.size() > 0
        
        # Tessellate!
        self.tessellate(tessellation)

    #-------------------------------------------------------------------------------
    # Radial tessellation utility function
    # Generate n-fold symmetric radial points, partition and store them, and 
    # tessellate
    #-------------------------------------------------------------------------------
    def computeRadialTessellation(self, tessellation, Nrings, 
                                  center             = None, 
                                  centralGenerator   = True,
                                  nfold              = 6,
                                  smallEdgeTolerance = 0.0):
        if center is None:
            center = [0.5*(self.xlimits[0] + self.xlimits[1]),
                      0.5*(self.ylimits[0] + self.ylimits[1])]
        assert (center[0] > self.xlimits[0] and center[0] < self.xlimits[1] and
                center[1] > self.ylimits[0] and center[1] < self.ylimits[1])
        
        # Radial domain sizing
        maxDistance = max(center[0] - self.xlimits[0], self.xlimits[1] - center[0],
                          center[1] - self.ylimits[0], self.ylimits[1] - center[1])
        dr = maxDistance/Nrings

        # Is the center point a generator?
        pps = []
        if centralGenerator:
            pps.append((center[0], center[1]))
            start   = 0.0
        else:
            start = 0.5

        # Create the points with N-fold symmetry
        for i in xrange(Nrings):
            r = (i + start)*dr
            nArcs = nfold*i
            for j in xrange(nArcs):
                theta = 2*pi*j/nArcs
                pps.append((center[0] + r*cos(theta), center[1] + r*sin(theta)))

        # Partition the generators
        self.points = self.partitioner.computePartition(pps)

        # Make sure boundary data is non-empty
        assert not self.PLC.empty()
        assert self.PLCpoints.size() > 0

        # Tessellate!
        self.tessellate(tessellation)

        if smallEdgeTolerance:
            self.cleanTessellation(tessellation, smallEdgeTolerance)

    #-------------------------------------------------------------------------------
    # Honeycomb tessellation utility function
    # Generate points that form regular hexagons, partition and store them, and
    # tessellate.
    #-------------------------------------------------------------------------------
    def computeHoneycombTessellation(self, tessellation, nx, ny):
        
        # Modifed point spacing
        dx = (self.xlimits[1] - self.xlimits[0])/(nx+1)
        dy = (self.ylimits[1] - self.ylimits[0])/(ny+1)
        xmin = self.xlimits[0]
        ymin = self.ylimits[0]

        pps = []
        iy = -1
        for ix in xrange(nx + 1):
            if (ix % 2) > 0:
                pps.append((xmin + (ix + 0.5)*dx,
                            ymin + (iy + 0.5 + 0.5*(ix % 2))*dy))

        for iy in xrange(ny + 1):
            for ix in xrange(nx + 1):
                pps.append((xmin + (ix + 0.5)*dx,
                            ymin + (iy + 0.5 + 0.5*(ix % 2))*dy))

        self.points = self.partitioner.computePartition(pps)

        # Make sure the boundary data corresponds to the specified limits
        boundaryPoints = boxListFromLimits(self.xlimits, self.ylimits)
        self.boundaryPoints = boundaryPoints
        self.PLC            = computePLCfromList2d(boundaryPoints)
        self.PLCpoints      = vectorize_double(boundaryPoints)

        # Tessellate!
        self.tessellate(tessellation)

    #-------------------------------------------------------------------------------
    # Generate a tessellation from a set of input points
    # Points are repartitioned (if asked), stored, and a mesh is output
    #-------------------------------------------------------------------------------
    def computeTessellationFromGenerators(self, tessellation, points, 
                                          boundaryPoints=None, 
                                          repartition=False):
        #assert len(points) != 0        
        if repartition:
            if ((type(points[0]) is not tuple) and (type(points[0]) is not list)):
                tupledPoints = []
                for i in range(len(points)/2):
                    tupledPoints.append((points[2*i], points[2*i+1]))
                points = tupledPoints
            self.points = self.partitioner.computePartition(points)
        else:
            self.points = vectorize_double(points)
                
        # Make sure the boundary data is non-empty. 
        # If a new boundary is specified, use that instead.
        if boundaryPoints is None:
            assert not self.PLC.empty()
            assert self.PLCpoints.size() > 0
        else:
            self.boundaryPoints = boundaryPoints
            self.PLC            = computePLCfromList2d(boundaryPoints)
            self.PLCpoints      = vectorize_double(boundaryPoints)

        # Tessellate!
        self.tessellate(tessellation)

    #-----------------------------------------------------------------------------
    # Apply N Lloyd iterations to the set of stored points to centroidally
    # smooth the output tessellation
    #-----------------------------------------------------------------------------
    def smoothTessellation(self, tessellation, iterations=10):
        smoothTessellation(tessellation,
                           self.tessellator,
                           self.points,
                           self.PLCpoints,
                           self.PLC,
                           iterations = iterations)

    #-----------------------------------------------------------------------------
    # Clean all edges of a tessellation below a relative tolerance per local zone
    #-----------------------------------------------------------------------------
    def cleanTessellation(self, tessellation, 
                          smallEdgeTolerance = 0.01,
                          meshEditor         = None):
        assert not tessellation.empty(), "Empty initial tessellation"
        assert smallEdgeTolerance > 0.0
        if meshEditor is None:
            meshEditor = polytope.MeshEditor2d(tessellation)
        meshEditor.cleanEdges(smallEdgeTolerance)    
    
    #-----------------------------------------------------------------------------
    # Output a silo file of a given tessellation having specified name and cycle
    # number (if part of a series). Generating points are stored by default.
    # Additional fields may be stored in the silo file by passing in dictionaries
    # of data on different centerings (i.e. cells, nodes, faces, etc.)
    #-----------------------------------------------------------------------------
    def outputTessellation(self, tessellation, name,
                           nodeFields = None,
                           edgeFields = None,
                           faceFields = None,
                           cellFields = None,
                           cycle      = 0):
        if cellFields is None:
            cellFields = {}
        genx = [self.points[2*i  ] for i in range(len(list(self.points))/2)]
        geny = [self.points[2*i+1] for i in range(len(list(self.points))/2)]
        cellFields["gen_x"] = genx
        cellFields["gen_y"] = geny
        polytope.writeTessellation2d(tessellation, name,
                                     nodeFields,
                                     edgeFields,
                                     faceFields,
                                     cellFields,
                                     cycle = cycle)
                   
