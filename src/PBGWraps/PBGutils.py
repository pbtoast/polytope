import sys
from pybindgen import *

#-------------------------------------------------------------------------------
# Add an object by name to a module/namespace, and publish it to the world.
#-------------------------------------------------------------------------------
def addObject(mod, name, *args, **kwargs):
    x = mod.add_class(name, *args, **kwargs)
    if not 'wrapObjs' in mod.__dict__:
        mod.wrapObjs = {}
    if "custom_template_class_name" in kwargs:
        pubname = kwargs["custom_template_class_name"]
    elif "python_name" in kwargs:
        pubname = kwargs["python_name"]
    else:
        pubname = name
    mod.wrapObjs[pubname] = x
    return x

#-------------------------------------------------------------------------------
# Add a reference symbol to a type.
#-------------------------------------------------------------------------------
def ref(name):
    return "%s&" % name

#-------------------------------------------------------------------------------
# Add a pointer symbol to a type.
#-------------------------------------------------------------------------------
def ptr(name):
    return "%s*" % name

def const_ptr(name):
    return "const %s*" % name

#-------------------------------------------------------------------------------
# Return the normal way we want to handle a pointer parameter.
#-------------------------------------------------------------------------------
def inputptrparam(cppobj, argname):
    return Parameter.new(ptr(cppobj), argname, transfer_ownership=False) # , direction=Parameter.DIRECTION_IN)

#-------------------------------------------------------------------------------
# Return the normal way we want to handle a reference parameter.
#-------------------------------------------------------------------------------
def refparam(cppobj, argname, default_value=None):
    return Parameter.new(ref(cppobj), argname, direction=Parameter.DIRECTION_INOUT, default_value=default_value)

def constrefparam(cppobj, argname, default_value=None):
    return Parameter.new(ref("const " + cppobj), argname, direction=Parameter.DIRECTION_INOUT, default_value=default_value)

#-------------------------------------------------------------------------------
# Generate the SWIG <-> pybindgen binding boilerplate.
# These methods can be used to write explicit SWIG in/out typemaps -- ugly
# but it does let pybindgen wrapped objects interact with SWIG wrapped objects.
#-------------------------------------------------------------------------------
def generateSWIGBindings(obj, out):
    pdict = {"class_name"   : obj.full_name,
             "pystruct"     : obj.get_pystruct(),
             "pytypestruct" : obj.pytypestruct}
    SWIGHelpersOut = '''

PyObject* pybindgen(%(class_name)s* obj) {
   %(pystruct)s *result;
   result = PyObject_New(%(pystruct)s, &%(pytypestruct)s);
   result->obj = obj;
   result->flags = PYBINDGEN_WRAPPER_FLAG_OBJECT_NOT_OWNED;
   return (PyObject*) result;
}

''' % pdict
    SWIGHelpersIn = '''

int unpybindgen(PyObject *o, %(class_name)s* obj) {
   POLY_ASSERT(o != NULL);
   POLY_ASSERT(obj != NULL);
   PyErr_Clear();
   if (o->ob_type != &%(pytypestruct)s) return 0;
   %(pystruct)s* py_obj = (%(pystruct)s*) o;
   if (py_obj->obj == NULL) return 0;
   obj = py_obj->obj;
   return 1;
}

''' % pdict

    out.writeln(SWIGHelpersOut)
    out.writeln(SWIGHelpersIn)
