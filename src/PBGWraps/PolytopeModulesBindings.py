# -*- mode: python -*-
#-------------------------------------------------------------------------------
# Master PolytopeModules file for generating the Polytope python bindings.
#-------------------------------------------------------------------------------
import os, sys
from pybindgen import *
from PBGutils import generateSWIGBindings

#-------------------------------------------------------------------------------
# This is the implementation of the multi-section factory that allows us to
# write multiple source files for compiling.
# 
# This is directly cribbed from Gustavo's example from the ns3 package at
# http://code.nsnam.org/ns-3-dev/file/e15adc7172f1/bindings/python/ns3modulegen.py
#-------------------------------------------------------------------------------
class PolytopeMultiSectionFactory(module.MultiSectionFactory):

    def __init__(self, main_file_name, modules):
        print "main file name: ", main_file_name
        print "modules:  ", modules
        self.basename, ext = os.path.splitext(main_file_name)
        super(PolytopeMultiSectionFactory, self).__init__()
        self.main_file_name = main_file_name
        self.main_sink = FileCodeSink(open(main_file_name, "wt"))
        self.header_name = self.basename + ".hh"
        header_file_name = os.path.join(os.path.dirname(self.main_file_name), 
                                        '.',
                                        self.header_name)
        self.header_sink = FileCodeSink(open(header_file_name, "wt"))
        self.section_sinks = {'__main__': self.main_sink}
        for module in modules:
            section_name = '%s_%s' % (self.basename, module.replace('-', '_'))
            file_name = os.path.join(os.path.dirname(self.main_file_name), "%s.C" % section_name)
            sink = FileCodeSink(open(file_name, "wt"))
            self.section_sinks[section_name] = sink

    def get_section_code_sink(self, section_name):
        return self.section_sinks[section_name]

    def get_main_code_sink(self):
        return self.main_sink

    def get_common_header_code_sink(self):
        return self.header_sink

    def get_common_header_include(self):
        return '"%s"' % self.header_name

    def close(self):
        self.header_sink.file.close()
        self.main_sink.file.close()
        for sink in self.section_sinks.itervalues():
            sink.file.close()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# A few useful global pybindgen settings.
settings.allow_subclassing = True
settings.deprecated_virtuals = False

# The set packages we're going to process.
pkgs_string = "CXXContainers MeshEditor PLC Tessellation Tessellator Functions"
pkgs = pkgs_string.split()

# Extract the desired output files.
outfile = "PolytopeModules.C"

#outfile = sys.argv[-1].replace("Bindings", "").split(".")[0]+".C"
print "Parsing pybindgen packages: ", pkgs

#-------------------------------------------------------------------------------
# Create the PolytopeModules module.
#-------------------------------------------------------------------------------
mod = Module("PolytopeModules")

# Strip the "HAVE_foo" statements off the input string and store as bools
assert len(sys.argv) == 5
mod.have_triangle      = (sys.argv[1].capitalize() == "True")
mod.have_tetgen        = (sys.argv[2].capitalize() == "True")
mod.have_boost_voronoi = (sys.argv[3].capitalize() == "True")
mod.have_mpi           = (sys.argv[4].capitalize() == "True")

# Go through each package and add its stuff to the module.
for p in pkgs:
    print "Generating types for %s" % p
    modname = "%sModule" % p
    exec("import %s" % modname)
    mod.begin_section("PolytopeModules_%s" % p)
    exec("%s = %s.%s(mod)" % (p, modname, p))
    mod.end_section("PolytopeModules_%s" % p)

# Now bind methods to the objects.
for p in pkgs:
    print "Binding methods for objects in %s" % p
    mod.begin_section("PolytopeModules_%s" % p)
    exec("%s.generateBindings(mod)" % p)
    mod.end_section("PolytopeModules_%s" % p)

#-------------------------------------------------------------------------------
# Generate the C code.
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    print "Generating C code."
    out = PolytopeMultiSectionFactory(outfile, pkgs)
    print "Section sinks:  ", out.section_sinks

    mod.generate(out, "PolytopeModules")

    for p in pkgs:
        sink = out.get_section_code_sink("PolytopeModules_%s" % p)
        for obj in eval(p).objs:
            generateSWIGBindings(obj, sink)

    out.close()
