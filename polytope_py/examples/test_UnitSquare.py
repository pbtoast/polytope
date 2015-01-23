from polytope_test_utilities import *
from numpy import *

def getPowerOfTwo(x):
    if x == 1: return 0
    n = 0
    result = 0
    while result != 1:
        n += 1
        result = x >> n
    return n

# =================================================================================

# Input parameters
partition = "Voronoi"  # Voronoi, Lattice
n         = 1000       # compute nxn lattice tessellation
nmin      = 50         # minimum n for nxn lattice
nmax      = 100        # maximum n for nxn lattice
testrange = True       # if True, Test from nmin to nmax; if False, test nxn only

# =================================================================================

# MPI stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Checks
assert partition in ("Voronoi", "Lattice")

# Initialize the partioner
if partition == "Lattice":
    comm = MPI.COMM_WORLD
    nprocs = comm.size
    assert nprocs & (nprocs-1) == 0, "Number of processors is not a power of 2!"
    n = getPowerOfTwo(nprocs)
    if n % 2 == 0:
        nprocx = 1 << n/2
        nprocy = 1 << n/2
    else:
        nprocx = 1 << (n+1)/2
        nprocy = 1 << (n-1)/2
    partitioner = LatticePartitioner(nprocx = nprocx,
                                     nprocy = nprocy)
else:
    partitioner = VoronoiPartitioner(lloydIterations = 50)


mesh = polytope.Tessellation2d()
T = PyTessellator(serialTessellator = "Boost",
                  xlimits           = [0.0, 1.0],
                  ylimits           = [0.0, 1.0],
                  partitioner       = partitioner,
                  timer             = True)

if testrange:
    for i,n in enumerate(arange(nmax-nmin+1) + nmin):
        if rank == 0:  print "Constructing %ix%i lattice mesh" % (n,n)
        T.computeLatticeTessellation(tessellation = mesh, nx=n, ny=n)
        T.outputTessellation(mesh, "UnitSquare_range", cycle=i)

else:
    T.computeLatticeTessellation(tessellation = mesh, nx=n, ny=n)
    T.outputTessellation(mesh, "UnitSquare_static")

