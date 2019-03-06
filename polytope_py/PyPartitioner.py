from numpy import *
from mpi4py import *
from PyPolytope_utilities import *
from polytope import *

#=============================================================================
# PyPartitioner
# Parent class for partitioning a group of points
#=============================================================================
class PyPartitioner:
    def __init__(self,
                 serialTessellator = None,
                 boundaryPoints    = None,
                 xlimits           = [0.0, 1.0],
                 ylimits           = [0.0, 1.0]):
        # Store the input tessellator
        if serialTessellator is None:
            self.serialTessellator = polytope.TriangleTessellator2d()
        else:
            self.serialTessellator = serialTessellator

        # Store the domain limits
        self.xlimits = xlimits
        self.ylimits = ylimits
        
        # Store the input boundary as vector_of_double PLCpoints and a faceted PLC
        # Default is to initialize x- and y-limits
        if boundaryPoints is None:
            boundaryPoints = boxListFromLimits(xlimits, ylimits)
        self.boundaryPoints = boundaryPoints
        
        # MPI bidness
        self.comm = MPI.COMM_WORLD
        return

    def computePartition(self, inputPoints):
        assert False, "No implementation"
        return

#=============================================================================
# VoronoiPartitioner
# Generate a Voronoi mesh from a collection of random seed locations,
# centroidally smooth them, and use the resulting graph as the partition
#=============================================================================
class VoronoiPartitioner(PyPartitioner):
    def __init__(self,
                 serialTessellator = None,
                 boundaryPoints    = None,
                 xlimits           = [0.0, 1.0],
                 ylimits           = [0.0, 1.0],
                 lloydIterations   = 100,
                 seed              = 102309823):
        self.iterations = lloydIterations
        random.seed(seed)
        PyPartitioner.__init__(self, 
                               serialTessellator, 
                               boundaryPoints, 
                               xlimits,
                               ylimits)
        return

    def computePartition(self, inputPoints):
        procPoints = self.generateSeedPoints()        
        points = vector_of_double()
        for i,gen in enumerate(inputPoints):
            owner = 0
            if self.comm.size > 1:
                minDist = distance2(gen, procPoints[0])
                for j,pt in enumerate(procPoints[1:]):
                    dist = distance2(gen, pt)
                    if dist < minDist:
                        minDist = dist
                        owner = j+1
            if owner == self.comm.rank:
                for g in gen:
                    points.push_back(g)
        assert points.size() > 2, "Partitioner assigned %i points to %i" % (points.size(), self.comm.rank)
        return points

    def generateSeedPoints(self):
        xmin = self.xlimits[0]
        xmax = self.xlimits[1]
        ymin = self.ylimits[0]
        ymax = self.ylimits[1]
        xmid = 0.5*(xmin + xmax)
        ymid = 0.5*(ymin + ymax)
        Lx   = (xmax-xmin)
        Ly   = (ymax-ymin)
        procPoints = []
        if self.comm.size > 1:
            # Boundary information to compute tessellation
            #boundaryPoints = boxList(xmin,xmax,ymin,ymax)
            PLCpoints      = vectorize_double(self.boundaryPoints)
            PLC            = computePLCfromList2d(self.boundaryPoints)

            # Random generator points
            seedGenerators = vector_of_double(2*self.comm.size)
            for i in xrange(self.comm.size):
                xseed = xmid + 0.1*Lx*(random.rand() - 0.5)
                yseed = ymid + 0.1*Ly*(random.rand() - 0.5)
                seedGenerators[2*i  ] = xseed
                seedGenerators[2*i+1] = yseed

            # Compute a tessellation
            tess = polytope.Tessellation2d()
            smoothTessellation(tessellation = tess,
                               tessellator  = self.serialTessellator,
                               points       = seedGenerators,
                               PLCpoints    = PLCpoints, 
                               PLC          = PLC,
                               iterations   = self.iterations)

            # Return the relaxed seed point locations
            for i in range(self.comm.size):
                procPoints.append([seedGenerators[2*i], seedGenerators[2*i+1]])
        return procPoints

#=============================================================================
# LatticePartitioner
# Generate an NxM lattice as the domain partition
#=============================================================================
class LatticePartitioner(PyPartitioner):
    def __init__(self,
                 serialTessellator = None,
                 boundaryPoints    = None,
                 xlimits           = [0.0, 1.0],
                 ylimits           = [0.0, 1.0],
                 nprocx            = None,
                 nprocy            = None):
        self.npx = nprocx
        self.npy = nprocy
        PyPartitioner.__init__(self, 
                               serialTessellator, 
                               boundaryPoints, 
                               xlimits,
                               ylimits)
        assert self.comm.size == self.npx * self.npy

        # Assign indices on the processor lattice
        self.ipx = self.comm.rank % self.npx
        self.ipy = (self.comm.rank - self.ipx) / self.npx
        assert self.ipx < self.npx
        assert self.ipy < self.npy

        # Assign physical bounds for each proc
        dx = (xlimits[1]-xlimits[0])/self.npx
        dy = (ylimits[1]-ylimits[0])/self.npy
        self.xlow  = xlimits[0] + dx * self.ipx
        self.xhigh = xlimits[0] + dx * (self.ipx + 1)
        self.ylow  = ylimits[0] + dy * self.ipy
        self.yhigh = ylimits[0] + dy * (self.ipy + 1)
        return

    def computePartition(self, inputPoints):
        points = vector_of_double()
        xpoints = array([g[0] for g in inputPoints])
        ypoints = array([g[1] for g in inputPoints])
        inds = where((xpoints >= self.xlow) * (xpoints < self.xhigh) *
                     (ypoints >= self.ylow) * (ypoints < self.yhigh))[0]
        for x,y in zip(xpoints[inds], ypoints[inds]):
            points.push_back(x)
            points.push_back(y)
        assert points.size() > 2
        return points

#=============================================================================
# RandomPartitioner
# No partition: each point gets assigned at random!
#=============================================================================
class RandomPartitioner(PyPartitioner):
    def __init__(self,
                 serialTessellator = None,
                 boundaryPoints    = None,
                 xlimits           = [0.0, 1.0],
                 ylimits           = [0.0, 1.0]):
        PyPartitioner.__init__(self, 
                               serialTessellator, 
                               boundaryPoints, 
                               xlimits,
                               ylimits)
        return

    def computePartition(self, inputPoints):
        points = vector_of_double()
        for i,gen in enumerate(inputPoints):
            owner = int(random.ran() * self.comm.size)
            if owner == self.comm.rank:
                for g in gen:
                    points.push_back(g)
        assert points.size() > 2
        return points
        
