from itertools import chain, izip
from numpy import array, where, pi, sin, cos, sqrt, prod

class VelocityField:
    def __init__(self,
                 L      = 1.0,
                 v0     = 1.0,
                 origin = [0.0, 0.0]):
        self.L      = L
        self.v0     = v0
        self.origin = origin
        return
    def update(self, points, dt):
        velocities = self.velocity(points)
        assert len(velocities) == len(points)
        newPoints = [p + 0.5*dt*v for (p,v) in zip(points, velocities)]
        velocities = self.velocity(newPoints)
        return [p + dt*v for (p,v) in zip(points, velocities)]
    def velocity(self, points):
        return [0.0 for p in points]

class SolidRotation(VelocityField):
    def __init__(self,
                 L      = 1.0,
                 v0     = 1.0,
                 origin = [0.0, 0.0]):
        VelocityField.__init__(self, L, v0, origin)
        return
    def velocity(self, points):
        N = len(points)/2
        x = array([(points[2*i  ] - self.origin[0]) for i in range(N)])
        y = array([(points[2*i+1] - self.origin[1]) for i in range(N)])
        u = -self.v0 * y
        v =  self.v0 * x
        return list(chain.from_iterable(izip(u,v)))

class TaylorGreenVortex(VelocityField):
    def __init__(self,
                 L      = 1.0,
                 v0     = 1.0,
                 origin = [0.0, 0.0]):
        VelocityField.__init__(self, L, v0, origin)
        return
    def velocity(self, points):
        N = len(points)/2
        x = array([(points[2*i  ] - 0.0*self.origin[0])*(1.0*pi/self.L) for i in range(N)])
        y = array([(points[2*i+1] - 0.0*self.origin[1])*(1.0*pi/self.L) for i in range(N)])
        u =  self.v0 * sin(x) * cos(y)
        v = -self.v0 * cos(x) * sin(y)
        return list(chain.from_iterable(izip(u,v)))

class GreshoVortex(VelocityField):
    def __init__(self,
                 L      = 1.0,
                 v0     = 1.0,
                 origin = [0.0, 0.0]):
        VelocityField.__init__(self, L, v0, origin)
        return
    def greshoVel(self, r):
        assert prod(r >= 0.0)
        vr = 0.0*r
        vr = where((r > 0.0) * (r <= 0.2), 5.0        , vr)
        vr = where((r > 0.2) * (r <= 0.4), (2.0/r-5.0), vr)
        return vr
    def velocity(self, points):
        N  = len(points)/2
        r  = array([sqrt((points[2*i]-self.origin[0])**2 + (points[2*i+1]-self.origin[1])**2) for i in range(N)])
        vr = self.greshoVel(r)
        x = array([(points[2*i  ] - self.origin[0]) for i in range(N)])
        y = array([(points[2*i+1] - self.origin[1]) for i in range(N)])
        u = -self.v0 * vr * y
        v =  self.v0 * vr * x
        return list(chain.from_iterable(izip(u,v)))

class DeformationVortex(VelocityField):
    def __init__(self,
                 L      = 1.0,
                 v0     = 1.0,
                 origin = [0.0, 0.0]):
        VelocityField.__init__(self, L, v0, origin)
        return
    def velocity(self, points):
        N = len(points)/2
        x = array([(points[2*i  ] - self.origin[0])*(4.0*pi/self.L) for i in range(N)])
        y = array([(points[2*i+1] - self.origin[1])*(4.0*pi/self.L) for i in range(N)])
        u =  self.v0 * sin(x) * cos(y)
        v = -self.v0 * cos(x) * sin(y)
        return list(chain.from_iterable(izip(u,v)))
