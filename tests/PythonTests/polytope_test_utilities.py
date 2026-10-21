import polytope
import random, builtins
from time import perf_counter
from contextlib import contextmanager
import Boundary2d

# Print statement wrapper for printing only on the root rank
def rootprint(string):
    comm = polytope.Communicator.instance()
    if (comm.getRank() == comm.getRoot()):
        print(string)

def outputMesh2d(mesh,
                 filePrefix,
                 fields=None,
                 cycle=0,
                 time=0.,
                 numFiles=-1):
    writer = polytope.SiloWriter2d(mesh)
    writer.generateTestVars()
    if (fields):
        for key, vals in fields.item():
            writer.addField(polytope.FieldCentering.Cell,
                            key,
                            vals)
    writer.write(filePrefix, cycle, time, numFiles)

class PointGenerator:
    def __init__(self, dim = 2,
                 boundary2d = None,
                 lobounds = None,
                 hibounds = None):
        self.dim = dim
        self.boundary2d = boundary2d
        xmin = None
        xmax = None
        if boundary2d:
            lobounds = boundary2d.low
            hibounds = boundary2d.high
        if lobounds and hibounds:
            if dim == 3:
                xmin = polytope.Point3d(*lobounds)
                xmax = polytope.Point3d(*hibounds)
            else:
                xmin = polytope.Point2d(*lobounds)
                xmax = polytope.Point2d(*hibounds)
        if (not xmax or not xmin):
            raise ValueError("Must specify lobounds/hibounds, or a boundary2d")
        self.xmin = xmin
        self.xmax = xmax

    def generate_random_points(self, N, seed = -1):
        if (seed >= 0):
            random.seed(seed)
        L = (self.xmax - self.xmin)
        pout = []
        if (self.boundary2d and self.boundary2d.mType > 1):
            while(len(pout) < self.dim*N):
                point = [self.xmin[d] + random.random()*L[d] for d in range(self.dim)]
                if (boundary2d.testInside(point)):
                    pout.extend(point)
        else:
            for _ in range(N):
                pout.extend([self.xmin[d] + random.random()*L[d] for d in range(self.dim)])
        return pout

    def generate_normal_random_points(self, N, seed = -1,
                                      mu = 0.5, sigma = 1./6.17737):
        if (seed >= 0):
            random.seed(seed)
        L = (self.xmax - self.xmin)
        pout = []
        if (self.boundary2d and self.boundary2d.mType > 1):
            while(len(pout) < self.dim*N):
                point = [self.xmin[d] + random.normalvariate(mu, sigma)*L[d] for d in range(self.dim)]
                if (boundary2d.testInside(point)):
                    pout.extend(point)
        else:
            while(len(pout) < self.dim*N):
                rnums = [random.normalvariate(mu, sigma) for _ in range(self.dim)]
                if (min(rnums) >= 0. and max(rnums) <= 1.):
                    for d, rnum in enumerate(rnums):
                        pout.append(self.xmin[d] + rnum*L[d])
        return pout

@contextmanager
def timer(name):
    """
    Time a code section. Run as:

    with timer("my_functions"):
        function_one()
        function_two()
    """
    comm = polytope.Communicator.instance()
    rank = comm.getRank()
    root = comm.getRoot()
    comm.Barrier()
    tstart = perf_counter()

    try:
        yield
    finally:
        comm.Barrier()
        elapsed = perf_counter() - tstart
        rootprint(f"{name}: {elapsed:.6f} seconds")
