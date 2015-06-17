from pybindgen import *
from PBGutils import *

#-------------------------------------------------------------------------------
# Helper method to add std::vector.
#-------------------------------------------------------------------------------
def generateStdVectorBindings(v, value, cppname, 
                              indexAsPointer = False,
                              wrapIterators = False):
    pointerValue = (value[-1:] == "*")

    # Constructors.
    v.add_constructor([])
    v.add_constructor([param("int", "size")])
    v.add_constructor([constrefparam(cppname, "rhs")])
    if not pointerValue:
        v.add_constructor([param("int", "size"), param(value, "value")])

    # __len__
    v.add_method("size", "unsigned int", [])
    v.add_method("size", "unsigned int", [], custom_name = "__len__")
    v.add_method("resize", None, [param("int", "size")])

    # __add__ and __iadd__
    v.add_function_as_method("concatContainers", cppname,
                             [param(cppname, "self"), param(cppname, "rhs")],
                             template_parameters = [cppname],
                             foreign_cpp_namespace = "polytope",
                             custom_name = "__add__")
    v.add_function_as_method("concatContainersInPlace", cppname,
                             [param(cppname, "self"), param(cppname, "rhs")],
                             template_parameters = [cppname],
                             foreign_cpp_namespace = "polytope",
                             custom_name = "__iadd__")

    # __mul__ and __imul__
    v.add_function_as_method("repeatContainer", cppname,
                             [param(cppname, "self"), param("unsigned int", "count")],
                             template_parameters = [cppname],
                             foreign_cpp_namespace = "polytope",
                             custom_name = "__mul__")
    v.add_function_as_method("repeatContainerInPlace", cppname,
                             [param(cppname, "self"), param("unsigned int", "count")],
                             template_parameters = [cppname],
                             foreign_cpp_namespace = "polytope",
                             custom_name = "__imul__")

    # __getitem__
    if pointerValue:
        v.add_function_as_method("indexContainer",
                                 retval(value, reference_existing_object=True),
                                 [param(cppname, "self"), param("int", "index")],
                                 template_parameters = [cppname],
                                 foreign_cpp_namespace = "polytope",
                                 custom_name = "__getitem__")
    else:
        if indexAsPointer:
            v.add_function_as_method("indexContainerAsPointer",
                                     retval(ptr(value), reference_existing_object=True),
                                     [param(cppname, "self"), param("int", "index")],
                                     template_parameters = [cppname],
                                     foreign_cpp_namespace = "polytope",
                                     custom_name = "__getitem__")
        else:
            v.add_function_as_method("indexContainer", value,
                                     [param(cppname, "self"), param("int", "index")],
                                     template_parameters = [cppname],
                                     foreign_cpp_namespace = "polytope",
                                     custom_name = "__getitem__")

    # __setitem__
    if not pointerValue:
        v.add_function_as_method("assignToContainerIndex", "int",
                                 [param(cppname, "self"), param("int", "index"), param(value, "value")],
                                 template_parameters = [cppname],
                                 foreign_cpp_namespace = "polytope",
                                 custom_name = "__setitem__")

    # append
    if pointerValue:
        v.add_function_as_method("appendToContainerOfPointers", "int",
                                 [param(cppname, "self"), param(value, "value", transfer_ownership=False)],
                                 template_parameters = [cppname],
                                 foreign_cpp_namespace = "polytope",
                                 custom_name = "append")
    else:
        v.add_method("push_back", None, [param(value, "value")])
        v.add_method("push_back", None, [param(value, "value")], custom_name="append")


    # __getslice__ and __setslice__.
    v.add_function_as_method("sliceContainer", cppname,
                             [param(cppname, "self"), param("int", "index1"), param("int", "index2")],
                             template_parameters = [cppname],
                             foreign_cpp_namespace = "polytope",
                             custom_name = "__getslice__")
    v.add_function_as_method("assignToSlice", "int",
                             [param(cppname, "self"), param("int", "index1"), param("int", "index2"), refparam(cppname, "values")],
                             template_parameters = [cppname],
                             foreign_cpp_namespace = "polytope",
                             custom_name = "__setslice__")
                                      
    # __contains__
    if pointerValue:
        v.add_function_as_method("containsPtr", "int",
                                 [param(cppname, "self"), param(value, "value", transfer_ownership=False)],
                                 template_parameters = [cppname],
                                 foreign_cpp_namespace = "polytope",
                                 custom_name = "__contains__")
    else:
        v.add_function_as_method("containsValue", "int",
                                 [param(cppname, "self"), param(value, "value")],
                                 template_parameters = [cppname],
                                 foreign_cpp_namespace = "polytope",
                                 custom_name = "__contains__")

    return

#-------------------------------------------------------------------------------
# Helper method to add std::set.
#-------------------------------------------------------------------------------
def generateStdSetBindings(v, value, cppname, 
                           indexAsPointer = False,
                           wrapIterators = False):
    pointerValue = (value[-1:] == "*")

    # Constructors.
    v.add_constructor([])
    v.add_constructor([constrefparam(cppname, "rhs")])

    # __len__
    v.add_method("size", "unsigned int", [])
    v.add_method("size", "unsigned int", [], custom_name = "__len__")

    # __contains__
    if pointerValue:
        v.add_function_as_method("containsPtr", "int",
                                 [param(cppname, "self"), param(value, "value", transfer_ownership=False)],
                                 template_parameters = [cppname],
                                 foreign_cpp_namespace = "polytope",
                                 custom_name = "__contains__")
    else:
        v.add_function_as_method("containsValue", "int",
                                 [param(cppname, "self"), param(value, "value")],
                                 template_parameters = [cppname],
                                 foreign_cpp_namespace = "polytope",
                                 custom_name = "__contains__")

    return

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class CXXContainers:

    #---------------------------------------------------------------------------
    # Add all out stuff.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes
        mod.add_include('"Polytope_CXXTypes.hh"')

        # Namespace.
        std = mod.add_cpp_namespace("std")

        self.objs = []

        # These are the basic types we form vectors of.
        # We break this set into two types:
        #   elementTypes0 : things we want to index by value
        #   elementTypes1 : things we want to index by pointer
        self.elementTypes0 = ["unsigned", 
                              "int", 
                              "double"]
        self.elementTypes1 = ["vector_of_unsigned", 
                              "vector_of_int", 
                              "vector_of_double", 
                              "vector_of_vector_of_int", 
                              "set_of_unsigned"]

        # std::vector types.
        for name in (self.elementTypes0 + self.elementTypes1):
            exec('''
self.vector_of_%(name)s = addObject(mod, "vector_of_%(name)s", allow_subclassing=True)
self.objs.append(self.vector_of_%(name)s)
''' % {"name" : name})

        # std::set types
        for name in self.elementTypes0:
            exec('''
self.set_of_%(name)s = addObject(mod, "set_of_%(name)s", allow_subclassing=True)
self.objs.append(self.set_of_%(name)s)
''' % {"name" : name})

        # A few value types need special mangling for the element specs.
        self.valueMap = {}
        for x in (self.elementTypes0 + self.elementTypes1):
            exec('self.valueMap["%s"] = "%s"' % (x, x))
        self.valueMap["unsigned"] = "unsigned int"
        self.valueMap["string"] = "std::string"
        self.valueMap["ULL"] = "uint64_t"

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        # Vector types.
        for name in self.elementTypes0:
            exec('generateStdVectorBindings(self.vector_of_%s, "%s", "vector_of_%s")' % (name, self.valueMap[name], name))
        for name in self.elementTypes1:
            exec('generateStdVectorBindings(self.vector_of_%s, "%s", "vector_of_%s", indexAsPointer=True)' % (name, name, name))

        # Set types.
        for name in self.elementTypes0:
            exec('generateStdSetBindings(self.set_of_%s, "%s", "set_of_%s")' % (name, self.valueMap[name], name))

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["std"]

