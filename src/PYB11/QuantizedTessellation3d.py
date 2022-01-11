from PYB11Generator import *

@PYB11template("IntType", "RealType")
class QuantizedTessellation2d:
    """QuantizedTessellation2d

An intermediate representation for 2D tessellations in integer
coordinates."""

    PYB11typedefs = """
  typedef Point3<%(IntType)s> IntPoint;
  typedef Point3<%(RealType)s> RealPoint;
"""

    #...........................................................................
    # Constructors
    def pyinit1(self,
                points = "const std::vector<%(RealType)s>&",
                boundaryPoints = "const std::vector<%(RealType)s>&"):
        "Construct with the given generators.  Finds the bounding limits and sets the quantized generators."

    @PYB11implementation("""[](const std::vector<%(RealType)s>& points,
                               py::tuple xmin_in,
                               py::tuple xmax_in) {
                                   %(RealType)s xmin[3] = {xmin_in[0].cast<%(RealType)s>(),
                                                           xmin_in[1].cast<%(RealType)s>(),
                                                           xmin_in[2].cast<%(RealType)s>()};
                                   %(RealType)s xmax[3] = {xmax_in[0].cast<%(RealType)s>(),
                                                           xmax_in[1].cast<%(RealType)s>(),
                                                           xmax_in[2].cast<%(RealType)s>()};
                                   return new QuantizedTessellation3d<%(IntType)s, %(RealType)s>(points, xmin, xmax);
                               }""")
    def pyinit2(self,
                points = "const std::vector<%(RealType)s>&",
                xmin = "py::tuple",
                xmax = "py::tuple"):
        "Construct with the given generators using the specified bounds."

    #...........................................................................
    # Methods
    @PYB11const
    @PYB11implementation("""[](const QuantizedTessellation3d<%(IntType)s, %(RealType)s>& self,
                               py::tuple& realcoords) {
                                   %(RealType)s rcoords[3] = {realcoords[0].cast<%(RealType)s>(),
                                                              realcoords[1].cast<%(RealType)s>(),
                                                              realcoords[2].cast<%(RealType)s>()};
                                   %(IntType)s icoords[3];
                                   self.quantize(rcoords, icoords);
                                   return py::make_tuple(icoords[0], icoords[1], icoords[2]);
                               }""")
    def quantize(self,
                 realcoords = "py::tuple"):
        "Convert real coordinates to integers."
        return "py::tuple"

    @PYB11const
    @PYB11implementation("""[](const QuantizedTessellation3d<%(IntType)s, %(RealType)s>& self,
                               py::tuple& intcoords) {
                                   %(IntType)s icoords[3] = {intcoords[0].cast<%(IntType)s>(),
                                                             intcoords[1].cast<%(IntType)s>(),
                                                             intcoords[2].cast<%(IntType)s>()};
                                   %(RealType)s rcoords[3];
                                   self.dequantize(icoords, rcoords);
                                   return py::make_tuple(rcoords[0], rcoords[1], rcoords[2]);
                               }""")
    def dequantize(self,
                   intcoords = "py::tuple"):
        "Convert int coordinates to real."
        return "py::tuple"

    @PYB11const
    def fillTessellation(self,
                         mesh = "Tessellation<3, %(RealType)s>&"):
        "Read out the current QuantizedTessellation to regular Tessellation."
        return "void"

    #...........................................................................
    # Properties
    xmin = PYB11property("py::tuple",
                         getterraw = """[](const QuantizedTessellation3d<%(IntType)s, %(RealType)s>& self) {
                             return py::make_tuple(self.xmin[0], self.xmin[1]);
                           }""",
                         setterraw = """[](QuantizedTessellation3d<%(IntType)s, %(RealType)s>& self,
                                           py::tuple val) {
                             self.xmin[0] = val[0].cast<%(RealType)s>();
                             self.xmin[1] = val[1].cast<%(RealType)s>();
                           }""",
                         doc = "The minimum (real) coordinate")

    xmax = PYB11property("py::tuple",
                         getterraw = """[](const QuantizedTessellation3d<%(IntType)s, %(RealType)s>& self) {
                             return py::make_tuple(self.xmax[0], self.xmax[1]);
                           }""",
                         setterraw = """[](QuantizedTessellation3d<%(IntType)s, %(RealType)s>& self,
                                           py::tuple val) {
                             self.xmax[0] = val[0].cast<%(RealType)s>();
                             self.xmax[1] = val[1].cast<%(RealType)s>();
                           }""",
                         doc = "The maximum (real) coordinate")

    #...........................................................................
    # Attributes
    length = PYB11readwrite()
    infRadius = PYB11readwrite()
    generators = PYB11readwrite()
    guardGenerators = PYB11readwrite()
    nodes = PYB11readwrite()
    edges = PYB11readwrite()
    faces = PYB11readwrite()
    cellFaces = PYB11readwrite()

#-------------------------------------------------------------------------------
# Template instantiations
QuantizedTessellation3d_KD = PYB11TemplateClass(QuantizedTessellation3d, template_parameters=("polytope::KeyTraits::Key", "double"), pyname="QuantizedTessellation3d")
