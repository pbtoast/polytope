"Master module for polytope"

from PYB11Generator import *

PYB11includes = ['"polytope.hh"',
                 '"polytope_write_OOGL.hh"',
                 '"polytope_pybind11_helpers.hh"']

PYB11namespaces = ["polytope"]

#-------------------------------------------------------------------------------
# Add the polytope classes
#-------------------------------------------------------------------------------
from PLC import *
from ReducedPLC import *
from Tessellation import *
from Tessellator import *
from Tessellators import *

#-------------------------------------------------------------------------------
# Polytope functions
#-------------------------------------------------------------------------------
@PYB11template("RealType")
def writePLCtoOFF(plc = "const PLC<%(Dimension)s, %(RealType)s>&",
                  coords = "const std::vector<%(RealType)s>&",
                  filename = "const std::string"):
    "Write a %(Dimension)sD OFF polylist file."
    return "void"

writePLCtoOFF2d = PYB11TemplateFunction(writePLCtoOFF, pyname="writePLCtoOFF", template_parameters={"RealType": "double",
                                                                                                    "Dimension": "2"})
writePLCtoOFF3d = PYB11TemplateFunction(writePLCtoOFF, pyname="writePLCtoOFF", template_parameters={"RealType": "double",
                                                                                                    "Dimension": "3"})

@PYB11template("int Dimension", "RealType")
@PYB11implementation("""[](const Tessellation<%(Dimension)s, %(RealType)s>& mesh,
                           std::string filePrefix,
                           py::dict nodeFieldsDict,
                           py::dict edgeFieldsDict,
                           py::dict faceFieldsDict,
                           py::dict cellFieldsDict,
                           int cycle,
                           %(RealType)s time) {
                               std::map<std::string, %(RealType)s*> nodeFields, edgeFields, faceFields, cellFields;
                               std::vector<std::vector<%(RealType)s>> nodeVals, edgeVals, faceVals, cellVals;
                               for (const auto item: nodeFieldsDict) {
                                 nodeVals.push_back(std::vector<double>());
                                 const auto key = item.first.cast<std::string>();
                                 py::list vals = item.second.cast<py::list>();
                                 for (const auto vali: vals) nodeVals.back().push_back(vali.cast<double>());
                                 nodeFields[key] = &nodeVals.back().front();
                               }
                               SiloWriter<%(Dimension)s, %(RealType)s>::write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, cycle, time);
                           }""")
def writeTessellation(mesh = "const Tessellation<%(Dimension)s, %(RealType)s>&",
                      filePrefix = "std::string",
                      nodeFieldsDict = "py::dict",
                      edgeFieldsDict = "py::dict",
                      faceFieldsDict = "py::dict",
                      cellFieldsDict = "py::dict",
                      cycle = "int",
                      time = "%(RealType)s"):
    "Write a tessellation (and arbitrary fields) to a silo file."
    return "void"

writeTessellation2d = PYB11TemplateFunction(writeTessellation, template_parameters=("2", "double"))
