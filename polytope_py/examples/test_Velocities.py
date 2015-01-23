from polytope_test_utilities import *
from numpy import *
from VelocityFields import *

# =================================================================================

# Input params
nx  = 100
ny  = 100
xb  = [0.0, 1.0]
yb  = [0.0, 1.0]
v0  = 1.0
cfl = 0.5
cycles = 200
startCycle = 1
dumpEvery = 10
runTests = [1,2,3,4]     # 1: Solid Rotation, 2: Taylor-Green, 3: Gresho, 4: 4x4 Vortices
serialTess = "Triangle"

# Derived params
L         = xb[1] - xb[0]
assert L == (yb[1] - yb[0])
origin    = [0.5*(xb[0]+xb[1]), 0.5*(yb[0]+yb[1])]
dx        = L/nx
dy        = L/ny
maxDT     = min(dx,dy)/v0
DT        = cfl*maxDT
testNames = ["SolidRotation", "TaylorGreenVortex", "GreshoVortex", "DeformationVortex"]

# MPI stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

assert serialTess in ("Boost", "Triangle")

mesh = polytope.Tessellation2d()
T = PyTessellator(serialTessellator = serialTess
                  xlimits = xb,
                  ylimits = yb)

# =================================================================================

for i, itest in enumerate(runTests):
    name = testNames[itest-1]
    outputName = "Velocities_%s" % name
    exec("V = %s(L = L, v0 = v0, origin = origin)" % name)
    T.computeLatticeTessellation(tessellation = mesh, nx=nx, ny=ny)
    T.outputTessellation(mesh, outputName, cycle=0)
    boundary = [xb[0], yb[0], xb[1], yb[0], xb[1], yb[1], xb[0], yb[1]]
    if rank == 0:  print "\nTest %i - %s" % (itest, name)
    for j in range(cycles):
        icycle = j+1
        if rank == 0:  print "  Cycle %i of %i" % (icycle, cycles)
        newPoints = V.update(T.points, DT)
        if itest == 1:  boundary  = V.update(boundary, DT)
        if icycle >= startCycle:
            T.computeTessellationFromGenerators(mesh, newPoints, boundaryPoints=boundary)
            if (icycle % dumpEvery == 0): T.outputTessellation(mesh, outputName, cycle=icycle)
        else:
            vecNewPoints = vector_of_double(len(newPoints))
            for i,p in enumerate(newPoints):  vecNewPoints[i] = p
            T.points = vecNewPoints
        
