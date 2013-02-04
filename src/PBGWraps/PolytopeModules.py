#include "PolytopeModules.hh"
static PyMethodDef PolytopeModules_polytope_functions[] = {
    {NULL, NULL, 0, NULL}
};
/* --- classes --- */


static PyObject *
initPolytopeModules_polytope(void)
{
    PyObject *m;
    m = Py_InitModule3((char *) "PolytopeModules.polytope", PolytopeModules_polytope_functions, NULL);
    if (m == NULL) {
        return NULL;
    }
    /* Register the 'polytope::PLC2d' class */
    if (PyType_Ready(&PyPolytopePLC2d_Type)) {
        return NULL;
    }
    PyModule_AddObject(m, (char *) "PLC2d", (PyObject *) &PyPolytopePLC2d_Type);
    /* Register the 'polytope::PLC3d' class */
    if (PyType_Ready(&PyPolytopePLC3d_Type)) {
        return NULL;
    }
    PyModule_AddObject(m, (char *) "PLC3d", (PyObject *) &PyPolytopePLC3d_Type);
    /* Register the 'polytope::Tessellation2d' class */
    if (PyType_Ready(&PyPolytopeTessellation2d_Type)) {
        return NULL;
    }
    PyModule_AddObject(m, (char *) "Tessellation2d", (PyObject *) &PyPolytopeTessellation2d_Type);
    /* Register the 'polytope::Tessellation3d' class */
    if (PyType_Ready(&PyPolytopeTessellation3d_Type)) {
        return NULL;
    }
    PyModule_AddObject(m, (char *) "Tessellation3d", (PyObject *) &PyPolytopeTessellation3d_Type);
    /* Register the 'polytope::Tessellator2d' class */
    if (PyType_Ready(&PyPolytopeTessellator2d_Type)) {
        return NULL;
    }
    PyModule_AddObject(m, (char *) "Tessellator2d", (PyObject *) &PyPolytopeTessellator2d_Type);
    /* Register the 'polytope::Tessellator3d' class */
    if (PyType_Ready(&PyPolytopeTessellator3d_Type)) {
        return NULL;
    }
    PyModule_AddObject(m, (char *) "Tessellator3d", (PyObject *) &PyPolytopeTessellator3d_Type);
    /* Register the 'polytope::TriangleTessellator' class */
    PyPolytopeTriangleTessellator_Type.tp_base = &PyPolytopeTessellator2d_Type;
    if (PyType_Ready(&PyPolytopeTriangleTessellator_Type)) {
        return NULL;
    }
    PyModule_AddObject(m, (char *) "TriangleTessellator", (PyObject *) &PyPolytopeTriangleTessellator_Type);
    /* Register the 'polytope::TetgenTessellator' class */
    PyPolytopeTetgenTessellator_Type.tp_base = &PyPolytopeTessellator3d_Type;
    if (PyType_Ready(&PyPolytopeTetgenTessellator_Type)) {
        return NULL;
    }
    PyModule_AddObject(m, (char *) "TetgenTessellator", (PyObject *) &PyPolytopeTetgenTessellator_Type);
    return m;
}
static PyMethodDef PolytopeModules_functions[] = {
    {NULL, NULL, 0, NULL}
};
/* --- containers --- */



static void
Pystd__set__lt__unsigned__gt__Iter__tp_clear(Pystd__set__lt__unsigned__gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

}


static int
Pystd__set__lt__unsigned__gt__Iter__tp_traverse(Pystd__set__lt__unsigned__gt__Iter *self, visitproc visit, void *arg)
{
    Py_VISIT((PyObject *) self->container);
    return 0;
}


static void
_wrap_Pystd__set__lt__unsigned__gt____tp_dealloc(Pystd__set__lt__unsigned__gt__ *self)
{
    delete self->obj;
    self->obj = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static void
_wrap_Pystd__set__lt__unsigned__gt__Iter__tp_dealloc(Pystd__set__lt__unsigned__gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static PyObject*
_wrap_Pystd__set__lt__unsigned__gt____tp_iter(Pystd__set__lt__unsigned__gt__ *self)
{
    Pystd__set__lt__unsigned__gt__Iter *iter = PyObject_GC_New(Pystd__set__lt__unsigned__gt__Iter, &Pystd__set__lt__unsigned__gt__Iter_Type);
    Py_INCREF(self);
    iter->container = self;
    iter->iterator = new std::set<unsigned>::iterator(self->obj->begin());
    return (PyObject*) iter;
}


static PyObject*
_wrap_Pystd__set__lt__unsigned__gt__Iter__tp_iter(Pystd__set__lt__unsigned__gt__Iter *self)
{
    Py_INCREF(self);
    return (PyObject*) self;
}

static PyObject* _wrap_Pystd__set__lt__unsigned__gt__Iter__tp_iternext(Pystd__set__lt__unsigned__gt__Iter *self)
{
    PyObject *py_retval;
    std::set<unsigned>::iterator iter;
    
    iter = *self->iterator;
    if (iter == self->container->obj->end()) {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    ++(*self->iterator);
    py_retval = Py_BuildValue((char *) "N", PyLong_FromUnsignedLong((*iter)));
    return py_retval;
}

int _wrap_convert_py2c__unsigned_int(PyObject *value, unsigned int *address)
{
    PyObject *py_retval;
    
    py_retval = Py_BuildValue((char *) "(O)", value);
    if (!PyArg_ParseTuple(py_retval, (char *) "I", &*address)) {
        Py_DECREF(py_retval);
        return 0;
    }
    Py_DECREF(py_retval);
    return 1;
}


int _wrap_convert_py2c__std__set__lt___unsigned___gt__(PyObject *arg, std::set<unsigned> *container)
{
    if (PyObject_IsInstance(arg, (PyObject*) &Pystd__set__lt__unsigned__gt___Type)) {
        *container = *((Pystd__set__lt__unsigned__gt__*)arg)->obj;
    } else if (PyList_Check(arg)) {
        container->clear();
        Py_ssize_t size = PyList_Size(arg);
        for (Py_ssize_t i = 0; i < size; i++) {
            unsigned int item;
            if (!_wrap_convert_py2c__unsigned_int(PyList_GET_ITEM(arg, i), &item)) {
                return 0;
            }
            container->insert(item);
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "parameter must be None, a set_of_uints instance, or a list of unsigned int");
        return 0;
    }
    return 1;
}


static int
_wrap_Pystd__set__lt__unsigned__gt____tp_init(Pystd__set__lt__unsigned__gt__ *self, PyObject *args, PyObject *kwargs)
{
    const char *keywords[] = {"arg", NULL};
    PyObject *arg = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, (char *) "|O", (char **) keywords, &arg)) {
        return -1;
    }

    self->obj = new std::set<unsigned>;

    if (arg == NULL)
        return 0;

    if (!_wrap_convert_py2c__std__set__lt___unsigned___gt__(arg, self->obj)) {
        delete self->obj;
        self->obj = NULL;
        return -1;
    }
    return 0;
}

PyTypeObject Pystd__set__lt__unsigned__gt___Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.set_of_uints",            /* tp_name */
    sizeof(Pystd__set__lt__unsigned__gt__),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__set__lt__unsigned__gt____tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)NULL,     /* tp_traverse */
    (inquiry)NULL,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__set__lt__unsigned__gt____tp_iter,          /* tp_iter */
    (iternextfunc)NULL,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)_wrap_Pystd__set__lt__unsigned__gt____tp_init,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};

PyTypeObject Pystd__set__lt__unsigned__gt__Iter_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.set_of_uintsIter",            /* tp_name */
    sizeof(Pystd__set__lt__unsigned__gt__Iter),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__set__lt__unsigned__gt__Iter__tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_HAVE_GC,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)Pystd__set__lt__unsigned__gt__Iter__tp_traverse,     /* tp_traverse */
    (inquiry)Pystd__set__lt__unsigned__gt__Iter__tp_clear,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__set__lt__unsigned__gt__Iter__tp_iter,          /* tp_iter */
    (iternextfunc)_wrap_Pystd__set__lt__unsigned__gt__Iter__tp_iternext,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)NULL,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};




static void
Pystd__vector__lt__int__gt__Iter__tp_clear(Pystd__vector__lt__int__gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

}


static int
Pystd__vector__lt__int__gt__Iter__tp_traverse(Pystd__vector__lt__int__gt__Iter *self, visitproc visit, void *arg)
{
    Py_VISIT((PyObject *) self->container);
    return 0;
}


static void
_wrap_Pystd__vector__lt__int__gt____tp_dealloc(Pystd__vector__lt__int__gt__ *self)
{
    delete self->obj;
    self->obj = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static void
_wrap_Pystd__vector__lt__int__gt__Iter__tp_dealloc(Pystd__vector__lt__int__gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static PyObject*
_wrap_Pystd__vector__lt__int__gt____tp_iter(Pystd__vector__lt__int__gt__ *self)
{
    Pystd__vector__lt__int__gt__Iter *iter = PyObject_GC_New(Pystd__vector__lt__int__gt__Iter, &Pystd__vector__lt__int__gt__Iter_Type);
    Py_INCREF(self);
    iter->container = self;
    iter->iterator = new std::vector<int>::iterator(self->obj->begin());
    return (PyObject*) iter;
}


static PyObject*
_wrap_Pystd__vector__lt__int__gt__Iter__tp_iter(Pystd__vector__lt__int__gt__Iter *self)
{
    Py_INCREF(self);
    return (PyObject*) self;
}

static PyObject* _wrap_Pystd__vector__lt__int__gt__Iter__tp_iternext(Pystd__vector__lt__int__gt__Iter *self)
{
    PyObject *py_retval;
    std::vector<int>::iterator iter;
    
    iter = *self->iterator;
    if (iter == self->container->obj->end()) {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    ++(*self->iterator);
    py_retval = Py_BuildValue((char *) "i", (*iter));
    return py_retval;
}

int _wrap_convert_py2c__int(PyObject *value, int *address)
{
    PyObject *py_retval;
    
    py_retval = Py_BuildValue((char *) "(O)", value);
    if (!PyArg_ParseTuple(py_retval, (char *) "i", &*address)) {
        Py_DECREF(py_retval);
        return 0;
    }
    Py_DECREF(py_retval);
    return 1;
}


int _wrap_convert_py2c__std__vector__lt___int___gt__(PyObject *arg, std::vector<int> *container)
{
    if (PyObject_IsInstance(arg, (PyObject*) &Pystd__vector__lt__int__gt___Type)) {
        *container = *((Pystd__vector__lt__int__gt__*)arg)->obj;
    } else if (PyList_Check(arg)) {
        container->clear();
        Py_ssize_t size = PyList_Size(arg);
        for (Py_ssize_t i = 0; i < size; i++) {
            int item;
            if (!_wrap_convert_py2c__int(PyList_GET_ITEM(arg, i), &item)) {
                return 0;
            }
            container->push_back(item);
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "parameter must be None, a vector_of_ints instance, or a list of int");
        return 0;
    }
    return 1;
}


static int
_wrap_Pystd__vector__lt__int__gt____tp_init(Pystd__vector__lt__int__gt__ *self, PyObject *args, PyObject *kwargs)
{
    const char *keywords[] = {"arg", NULL};
    PyObject *arg = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, (char *) "|O", (char **) keywords, &arg)) {
        return -1;
    }

    self->obj = new std::vector<int>;

    if (arg == NULL)
        return 0;

    if (!_wrap_convert_py2c__std__vector__lt___int___gt__(arg, self->obj)) {
        delete self->obj;
        self->obj = NULL;
        return -1;
    }
    return 0;
}

PyTypeObject Pystd__vector__lt__int__gt___Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_ints",            /* tp_name */
    sizeof(Pystd__vector__lt__int__gt__),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__int__gt____tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)NULL,     /* tp_traverse */
    (inquiry)NULL,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__int__gt____tp_iter,          /* tp_iter */
    (iternextfunc)NULL,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)_wrap_Pystd__vector__lt__int__gt____tp_init,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};

PyTypeObject Pystd__vector__lt__int__gt__Iter_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_intsIter",            /* tp_name */
    sizeof(Pystd__vector__lt__int__gt__Iter),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__int__gt__Iter__tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_HAVE_GC,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)Pystd__vector__lt__int__gt__Iter__tp_traverse,     /* tp_traverse */
    (inquiry)Pystd__vector__lt__int__gt__Iter__tp_clear,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__int__gt__Iter__tp_iter,          /* tp_iter */
    (iternextfunc)_wrap_Pystd__vector__lt__int__gt__Iter__tp_iternext,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)NULL,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};




static void
Pystd__vector__lt__unsigned__gt__Iter__tp_clear(Pystd__vector__lt__unsigned__gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

}


static int
Pystd__vector__lt__unsigned__gt__Iter__tp_traverse(Pystd__vector__lt__unsigned__gt__Iter *self, visitproc visit, void *arg)
{
    Py_VISIT((PyObject *) self->container);
    return 0;
}


static void
_wrap_Pystd__vector__lt__unsigned__gt____tp_dealloc(Pystd__vector__lt__unsigned__gt__ *self)
{
    delete self->obj;
    self->obj = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static void
_wrap_Pystd__vector__lt__unsigned__gt__Iter__tp_dealloc(Pystd__vector__lt__unsigned__gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static PyObject*
_wrap_Pystd__vector__lt__unsigned__gt____tp_iter(Pystd__vector__lt__unsigned__gt__ *self)
{
    Pystd__vector__lt__unsigned__gt__Iter *iter = PyObject_GC_New(Pystd__vector__lt__unsigned__gt__Iter, &Pystd__vector__lt__unsigned__gt__Iter_Type);
    Py_INCREF(self);
    iter->container = self;
    iter->iterator = new std::vector<unsigned>::iterator(self->obj->begin());
    return (PyObject*) iter;
}


static PyObject*
_wrap_Pystd__vector__lt__unsigned__gt__Iter__tp_iter(Pystd__vector__lt__unsigned__gt__Iter *self)
{
    Py_INCREF(self);
    return (PyObject*) self;
}

static PyObject* _wrap_Pystd__vector__lt__unsigned__gt__Iter__tp_iternext(Pystd__vector__lt__unsigned__gt__Iter *self)
{
    PyObject *py_retval;
    std::vector<unsigned>::iterator iter;
    
    iter = *self->iterator;
    if (iter == self->container->obj->end()) {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    ++(*self->iterator);
    py_retval = Py_BuildValue((char *) "N", PyLong_FromUnsignedLong((*iter)));
    return py_retval;
}

int _wrap_convert_py2c__std__vector__lt___unsigned___gt__(PyObject *arg, std::vector<unsigned> *container)
{
    if (PyObject_IsInstance(arg, (PyObject*) &Pystd__vector__lt__unsigned__gt___Type)) {
        *container = *((Pystd__vector__lt__unsigned__gt__*)arg)->obj;
    } else if (PyList_Check(arg)) {
        container->clear();
        Py_ssize_t size = PyList_Size(arg);
        for (Py_ssize_t i = 0; i < size; i++) {
            unsigned int item;
            if (!_wrap_convert_py2c__unsigned_int(PyList_GET_ITEM(arg, i), &item)) {
                return 0;
            }
            container->push_back(item);
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "parameter must be None, a vector_of_uints instance, or a list of unsigned int");
        return 0;
    }
    return 1;
}


static int
_wrap_Pystd__vector__lt__unsigned__gt____tp_init(Pystd__vector__lt__unsigned__gt__ *self, PyObject *args, PyObject *kwargs)
{
    const char *keywords[] = {"arg", NULL};
    PyObject *arg = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, (char *) "|O", (char **) keywords, &arg)) {
        return -1;
    }

    self->obj = new std::vector<unsigned>;

    if (arg == NULL)
        return 0;

    if (!_wrap_convert_py2c__std__vector__lt___unsigned___gt__(arg, self->obj)) {
        delete self->obj;
        self->obj = NULL;
        return -1;
    }
    return 0;
}

PyTypeObject Pystd__vector__lt__unsigned__gt___Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_uints",            /* tp_name */
    sizeof(Pystd__vector__lt__unsigned__gt__),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__unsigned__gt____tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)NULL,     /* tp_traverse */
    (inquiry)NULL,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__unsigned__gt____tp_iter,          /* tp_iter */
    (iternextfunc)NULL,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)_wrap_Pystd__vector__lt__unsigned__gt____tp_init,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};

PyTypeObject Pystd__vector__lt__unsigned__gt__Iter_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_uintsIter",            /* tp_name */
    sizeof(Pystd__vector__lt__unsigned__gt__Iter),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__unsigned__gt__Iter__tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_HAVE_GC,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)Pystd__vector__lt__unsigned__gt__Iter__tp_traverse,     /* tp_traverse */
    (inquiry)Pystd__vector__lt__unsigned__gt__Iter__tp_clear,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__unsigned__gt__Iter__tp_iter,          /* tp_iter */
    (iternextfunc)_wrap_Pystd__vector__lt__unsigned__gt__Iter__tp_iternext,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)NULL,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};




static void
Pystd__vector__lt__double__gt__Iter__tp_clear(Pystd__vector__lt__double__gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

}


static int
Pystd__vector__lt__double__gt__Iter__tp_traverse(Pystd__vector__lt__double__gt__Iter *self, visitproc visit, void *arg)
{
    Py_VISIT((PyObject *) self->container);
    return 0;
}


static void
_wrap_Pystd__vector__lt__double__gt____tp_dealloc(Pystd__vector__lt__double__gt__ *self)
{
    delete self->obj;
    self->obj = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static void
_wrap_Pystd__vector__lt__double__gt__Iter__tp_dealloc(Pystd__vector__lt__double__gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static PyObject*
_wrap_Pystd__vector__lt__double__gt____tp_iter(Pystd__vector__lt__double__gt__ *self)
{
    Pystd__vector__lt__double__gt__Iter *iter = PyObject_GC_New(Pystd__vector__lt__double__gt__Iter, &Pystd__vector__lt__double__gt__Iter_Type);
    Py_INCREF(self);
    iter->container = self;
    iter->iterator = new std::vector<double>::iterator(self->obj->begin());
    return (PyObject*) iter;
}


static PyObject*
_wrap_Pystd__vector__lt__double__gt__Iter__tp_iter(Pystd__vector__lt__double__gt__Iter *self)
{
    Py_INCREF(self);
    return (PyObject*) self;
}

static PyObject* _wrap_Pystd__vector__lt__double__gt__Iter__tp_iternext(Pystd__vector__lt__double__gt__Iter *self)
{
    PyObject *py_retval;
    std::vector<double>::iterator iter;
    
    iter = *self->iterator;
    if (iter == self->container->obj->end()) {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    ++(*self->iterator);
    py_retval = Py_BuildValue((char *) "d", (*iter));
    return py_retval;
}

int _wrap_convert_py2c__double(PyObject *value, double *address)
{
    PyObject *py_retval;
    
    py_retval = Py_BuildValue((char *) "(O)", value);
    if (!PyArg_ParseTuple(py_retval, (char *) "d", &*address)) {
        Py_DECREF(py_retval);
        return 0;
    }
    Py_DECREF(py_retval);
    return 1;
}


int _wrap_convert_py2c__std__vector__lt___double___gt__(PyObject *arg, std::vector<double> *container)
{
    if (PyObject_IsInstance(arg, (PyObject*) &Pystd__vector__lt__double__gt___Type)) {
        *container = *((Pystd__vector__lt__double__gt__*)arg)->obj;
    } else if (PyList_Check(arg)) {
        container->clear();
        Py_ssize_t size = PyList_Size(arg);
        for (Py_ssize_t i = 0; i < size; i++) {
            double item;
            if (!_wrap_convert_py2c__double(PyList_GET_ITEM(arg, i), &item)) {
                return 0;
            }
            container->push_back(item);
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "parameter must be None, a vector_of_doubles instance, or a list of double");
        return 0;
    }
    return 1;
}


static int
_wrap_Pystd__vector__lt__double__gt____tp_init(Pystd__vector__lt__double__gt__ *self, PyObject *args, PyObject *kwargs)
{
    const char *keywords[] = {"arg", NULL};
    PyObject *arg = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, (char *) "|O", (char **) keywords, &arg)) {
        return -1;
    }

    self->obj = new std::vector<double>;

    if (arg == NULL)
        return 0;

    if (!_wrap_convert_py2c__std__vector__lt___double___gt__(arg, self->obj)) {
        delete self->obj;
        self->obj = NULL;
        return -1;
    }
    return 0;
}

PyTypeObject Pystd__vector__lt__double__gt___Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_doubles",            /* tp_name */
    sizeof(Pystd__vector__lt__double__gt__),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__double__gt____tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)NULL,     /* tp_traverse */
    (inquiry)NULL,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__double__gt____tp_iter,          /* tp_iter */
    (iternextfunc)NULL,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)_wrap_Pystd__vector__lt__double__gt____tp_init,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};

PyTypeObject Pystd__vector__lt__double__gt__Iter_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_doublesIter",            /* tp_name */
    sizeof(Pystd__vector__lt__double__gt__Iter),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__double__gt__Iter__tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_HAVE_GC,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)Pystd__vector__lt__double__gt__Iter__tp_traverse,     /* tp_traverse */
    (inquiry)Pystd__vector__lt__double__gt__Iter__tp_clear,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__double__gt__Iter__tp_iter,          /* tp_iter */
    (iternextfunc)_wrap_Pystd__vector__lt__double__gt__Iter__tp_iternext,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)NULL,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};




static void
Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter__tp_clear(Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

}


static int
Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter__tp_traverse(Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter *self, visitproc visit, void *arg)
{
    Py_VISIT((PyObject *) self->container);
    return 0;
}


static void
_wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt____tp_dealloc(Pystd__vector__lt__std__vector__lt__int__gt_____gt__ *self)
{
    delete self->obj;
    self->obj = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static void
_wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter__tp_dealloc(Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static PyObject*
_wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt____tp_iter(Pystd__vector__lt__std__vector__lt__int__gt_____gt__ *self)
{
    Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter *iter = PyObject_GC_New(Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter, &Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter_Type);
    Py_INCREF(self);
    iter->container = self;
    iter->iterator = new std::vector<std::vector<int> >::iterator(self->obj->begin());
    return (PyObject*) iter;
}


static PyObject*
_wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter__tp_iter(Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter *self)
{
    Py_INCREF(self);
    return (PyObject*) self;
}

static PyObject* _wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter__tp_iternext(Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter *self)
{
    PyObject *py_retval;
    std::vector<std::vector<int> >::iterator iter;
    Pystd__vector__lt__int__gt__ *py_std__vector__lt__int__gt__;
    
    iter = *self->iterator;
    if (iter == self->container->obj->end()) {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    ++(*self->iterator);
    py_std__vector__lt__int__gt__ = PyObject_New(Pystd__vector__lt__int__gt__, &Pystd__vector__lt__int__gt___Type);
    py_std__vector__lt__int__gt__->obj = new std::vector<int>((*iter));
    py_retval = Py_BuildValue((char *) "N", py_std__vector__lt__int__gt__);
    return py_retval;
}

int _wrap_convert_py2c__std__vector__lt___std__vector__lt___int___gt_____gt__(PyObject *arg, std::vector<std::vector<int> > *container)
{
    if (PyObject_IsInstance(arg, (PyObject*) &Pystd__vector__lt__std__vector__lt__int__gt_____gt___Type)) {
        *container = *((Pystd__vector__lt__std__vector__lt__int__gt_____gt__*)arg)->obj;
    } else if (PyList_Check(arg)) {
        container->clear();
        Py_ssize_t size = PyList_Size(arg);
        for (Py_ssize_t i = 0; i < size; i++) {
            std::vector< int > item;
            if (!_wrap_convert_py2c__std__vector__lt___int___gt__(PyList_GET_ITEM(arg, i), &item)) {
                return 0;
            }
            container->push_back(item);
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "parameter must be None, a vector_of_vector_of_ints instance, or a list of std::vector< int >");
        return 0;
    }
    return 1;
}


static int
_wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt____tp_init(Pystd__vector__lt__std__vector__lt__int__gt_____gt__ *self, PyObject *args, PyObject *kwargs)
{
    const char *keywords[] = {"arg", NULL};
    PyObject *arg = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, (char *) "|O", (char **) keywords, &arg)) {
        return -1;
    }

    self->obj = new std::vector<std::vector<int> >;

    if (arg == NULL)
        return 0;

    if (!_wrap_convert_py2c__std__vector__lt___std__vector__lt___int___gt_____gt__(arg, self->obj)) {
        delete self->obj;
        self->obj = NULL;
        return -1;
    }
    return 0;
}

PyTypeObject Pystd__vector__lt__std__vector__lt__int__gt_____gt___Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_vector_of_ints",            /* tp_name */
    sizeof(Pystd__vector__lt__std__vector__lt__int__gt_____gt__),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt____tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)NULL,     /* tp_traverse */
    (inquiry)NULL,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt____tp_iter,          /* tp_iter */
    (iternextfunc)NULL,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)_wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt____tp_init,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};

PyTypeObject Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_vector_of_intsIter",            /* tp_name */
    sizeof(Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter__tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_HAVE_GC,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter__tp_traverse,     /* tp_traverse */
    (inquiry)Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter__tp_clear,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter__tp_iter,          /* tp_iter */
    (iternextfunc)_wrap_Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter__tp_iternext,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)NULL,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};




static void
Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter__tp_clear(Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

}


static int
Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter__tp_traverse(Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter *self, visitproc visit, void *arg)
{
    Py_VISIT((PyObject *) self->container);
    return 0;
}


static void
_wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt____tp_dealloc(Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__ *self)
{
    delete self->obj;
    self->obj = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static void
_wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter__tp_dealloc(Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static PyObject*
_wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt____tp_iter(Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__ *self)
{
    Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter *iter = PyObject_GC_New(Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter, &Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter_Type);
    Py_INCREF(self);
    iter->container = self;
    iter->iterator = new std::vector<std::vector<unsigned> >::iterator(self->obj->begin());
    return (PyObject*) iter;
}


static PyObject*
_wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter__tp_iter(Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter *self)
{
    Py_INCREF(self);
    return (PyObject*) self;
}

static PyObject* _wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter__tp_iternext(Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter *self)
{
    PyObject *py_retval;
    std::vector<std::vector<unsigned> >::iterator iter;
    Pystd__vector__lt__unsigned__gt__ *py_std__vector__lt__unsigned__gt__;
    
    iter = *self->iterator;
    if (iter == self->container->obj->end()) {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    ++(*self->iterator);
    py_std__vector__lt__unsigned__gt__ = PyObject_New(Pystd__vector__lt__unsigned__gt__, &Pystd__vector__lt__unsigned__gt___Type);
    py_std__vector__lt__unsigned__gt__->obj = new std::vector<unsigned>((*iter));
    py_retval = Py_BuildValue((char *) "N", py_std__vector__lt__unsigned__gt__);
    return py_retval;
}

int _wrap_convert_py2c__std__vector__lt___std__vector__lt___unsigned___gt_____gt__(PyObject *arg, std::vector<std::vector<unsigned> > *container)
{
    if (PyObject_IsInstance(arg, (PyObject*) &Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt___Type)) {
        *container = *((Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__*)arg)->obj;
    } else if (PyList_Check(arg)) {
        container->clear();
        Py_ssize_t size = PyList_Size(arg);
        for (Py_ssize_t i = 0; i < size; i++) {
            std::vector< unsigned > item;
            if (!_wrap_convert_py2c__std__vector__lt___unsigned___gt__(PyList_GET_ITEM(arg, i), &item)) {
                return 0;
            }
            container->push_back(item);
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "parameter must be None, a vector_of_vector_of_uints instance, or a list of std::vector< unsigned >");
        return 0;
    }
    return 1;
}


static int
_wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt____tp_init(Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__ *self, PyObject *args, PyObject *kwargs)
{
    const char *keywords[] = {"arg", NULL};
    PyObject *arg = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, (char *) "|O", (char **) keywords, &arg)) {
        return -1;
    }

    self->obj = new std::vector<std::vector<unsigned> >;

    if (arg == NULL)
        return 0;

    if (!_wrap_convert_py2c__std__vector__lt___std__vector__lt___unsigned___gt_____gt__(arg, self->obj)) {
        delete self->obj;
        self->obj = NULL;
        return -1;
    }
    return 0;
}

PyTypeObject Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt___Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_vector_of_uints",            /* tp_name */
    sizeof(Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt____tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)NULL,     /* tp_traverse */
    (inquiry)NULL,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt____tp_iter,          /* tp_iter */
    (iternextfunc)NULL,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)_wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt____tp_init,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};

PyTypeObject Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_vector_of_uintsIter",            /* tp_name */
    sizeof(Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter__tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_HAVE_GC,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter__tp_traverse,     /* tp_traverse */
    (inquiry)Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter__tp_clear,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter__tp_iter,          /* tp_iter */
    (iternextfunc)_wrap_Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter__tp_iternext,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)NULL,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};




static void
Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter__tp_clear(Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

}


static int
Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter__tp_traverse(Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter *self, visitproc visit, void *arg)
{
    Py_VISIT((PyObject *) self->container);
    return 0;
}


static void
_wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt____tp_dealloc(Pystd__vector__lt__std__vector__lt__double__gt_____gt__ *self)
{
    delete self->obj;
    self->obj = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static void
_wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter__tp_dealloc(Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static PyObject*
_wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt____tp_iter(Pystd__vector__lt__std__vector__lt__double__gt_____gt__ *self)
{
    Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter *iter = PyObject_GC_New(Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter, &Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter_Type);
    Py_INCREF(self);
    iter->container = self;
    iter->iterator = new std::vector<std::vector<double> >::iterator(self->obj->begin());
    return (PyObject*) iter;
}


static PyObject*
_wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter__tp_iter(Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter *self)
{
    Py_INCREF(self);
    return (PyObject*) self;
}

static PyObject* _wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter__tp_iternext(Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter *self)
{
    PyObject *py_retval;
    std::vector<std::vector<double> >::iterator iter;
    Pystd__vector__lt__double__gt__ *py_std__vector__lt__double__gt__;
    
    iter = *self->iterator;
    if (iter == self->container->obj->end()) {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    ++(*self->iterator);
    py_std__vector__lt__double__gt__ = PyObject_New(Pystd__vector__lt__double__gt__, &Pystd__vector__lt__double__gt___Type);
    py_std__vector__lt__double__gt__->obj = new std::vector<double>((*iter));
    py_retval = Py_BuildValue((char *) "N", py_std__vector__lt__double__gt__);
    return py_retval;
}

int _wrap_convert_py2c__std__vector__lt___std__vector__lt___double___gt_____gt__(PyObject *arg, std::vector<std::vector<double> > *container)
{
    if (PyObject_IsInstance(arg, (PyObject*) &Pystd__vector__lt__std__vector__lt__double__gt_____gt___Type)) {
        *container = *((Pystd__vector__lt__std__vector__lt__double__gt_____gt__*)arg)->obj;
    } else if (PyList_Check(arg)) {
        container->clear();
        Py_ssize_t size = PyList_Size(arg);
        for (Py_ssize_t i = 0; i < size; i++) {
            std::vector< double > item;
            if (!_wrap_convert_py2c__std__vector__lt___double___gt__(PyList_GET_ITEM(arg, i), &item)) {
                return 0;
            }
            container->push_back(item);
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "parameter must be None, a vector_of_vector_of_doubles instance, or a list of std::vector< double >");
        return 0;
    }
    return 1;
}


static int
_wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt____tp_init(Pystd__vector__lt__std__vector__lt__double__gt_____gt__ *self, PyObject *args, PyObject *kwargs)
{
    const char *keywords[] = {"arg", NULL};
    PyObject *arg = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, (char *) "|O", (char **) keywords, &arg)) {
        return -1;
    }

    self->obj = new std::vector<std::vector<double> >;

    if (arg == NULL)
        return 0;

    if (!_wrap_convert_py2c__std__vector__lt___std__vector__lt___double___gt_____gt__(arg, self->obj)) {
        delete self->obj;
        self->obj = NULL;
        return -1;
    }
    return 0;
}

PyTypeObject Pystd__vector__lt__std__vector__lt__double__gt_____gt___Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_vector_of_doubles",            /* tp_name */
    sizeof(Pystd__vector__lt__std__vector__lt__double__gt_____gt__),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt____tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)NULL,     /* tp_traverse */
    (inquiry)NULL,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt____tp_iter,          /* tp_iter */
    (iternextfunc)NULL,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)_wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt____tp_init,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};

PyTypeObject Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_vector_of_doublesIter",            /* tp_name */
    sizeof(Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter__tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_HAVE_GC,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter__tp_traverse,     /* tp_traverse */
    (inquiry)Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter__tp_clear,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter__tp_iter,          /* tp_iter */
    (iternextfunc)_wrap_Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter__tp_iternext,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)NULL,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};




static void
Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter__tp_clear(Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

}


static int
Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter__tp_traverse(Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter *self, visitproc visit, void *arg)
{
    Py_VISIT((PyObject *) self->container);
    return 0;
}


static void
_wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt____tp_dealloc(Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__ *self)
{
    delete self->obj;
    self->obj = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static void
_wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter__tp_dealloc(Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static PyObject*
_wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt____tp_iter(Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__ *self)
{
    Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter *iter = PyObject_GC_New(Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter, &Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter_Type);
    Py_INCREF(self);
    iter->container = self;
    iter->iterator = new std::vector<std::set<unsigned> >::iterator(self->obj->begin());
    return (PyObject*) iter;
}


static PyObject*
_wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter__tp_iter(Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter *self)
{
    Py_INCREF(self);
    return (PyObject*) self;
}

static PyObject* _wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter__tp_iternext(Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter *self)
{
    PyObject *py_retval;
    std::vector<std::set<unsigned> >::iterator iter;
    Pystd__set__lt__unsigned__gt__ *py_std__set__lt__unsigned__gt__;
    
    iter = *self->iterator;
    if (iter == self->container->obj->end()) {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    ++(*self->iterator);
    py_std__set__lt__unsigned__gt__ = PyObject_New(Pystd__set__lt__unsigned__gt__, &Pystd__set__lt__unsigned__gt___Type);
    py_std__set__lt__unsigned__gt__->obj = new std::set<unsigned>((*iter));
    py_retval = Py_BuildValue((char *) "N", py_std__set__lt__unsigned__gt__);
    return py_retval;
}

int _wrap_convert_py2c__std__vector__lt___std__set__lt___unsigned___gt_____gt__(PyObject *arg, std::vector<std::set<unsigned> > *container)
{
    if (PyObject_IsInstance(arg, (PyObject*) &Pystd__vector__lt__std__set__lt__unsigned__gt_____gt___Type)) {
        *container = *((Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__*)arg)->obj;
    } else if (PyList_Check(arg)) {
        container->clear();
        Py_ssize_t size = PyList_Size(arg);
        for (Py_ssize_t i = 0; i < size; i++) {
            std::set< unsigned > item;
            if (!_wrap_convert_py2c__std__set__lt___unsigned___gt__(PyList_GET_ITEM(arg, i), &item)) {
                return 0;
            }
            container->push_back(item);
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "parameter must be None, a vector_of_set_of_unsigned instance, or a list of std::set< unsigned >");
        return 0;
    }
    return 1;
}


static int
_wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt____tp_init(Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__ *self, PyObject *args, PyObject *kwargs)
{
    const char *keywords[] = {"arg", NULL};
    PyObject *arg = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, (char *) "|O", (char **) keywords, &arg)) {
        return -1;
    }

    self->obj = new std::vector<std::set<unsigned> >;

    if (arg == NULL)
        return 0;

    if (!_wrap_convert_py2c__std__vector__lt___std__set__lt___unsigned___gt_____gt__(arg, self->obj)) {
        delete self->obj;
        self->obj = NULL;
        return -1;
    }
    return 0;
}

PyTypeObject Pystd__vector__lt__std__set__lt__unsigned__gt_____gt___Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_set_of_unsigned",            /* tp_name */
    sizeof(Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt____tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)NULL,     /* tp_traverse */
    (inquiry)NULL,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt____tp_iter,          /* tp_iter */
    (iternextfunc)NULL,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)_wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt____tp_init,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};

PyTypeObject Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_set_of_unsignedIter",            /* tp_name */
    sizeof(Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter__tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_HAVE_GC,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter__tp_traverse,     /* tp_traverse */
    (inquiry)Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter__tp_clear,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter__tp_iter,          /* tp_iter */
    (iternextfunc)_wrap_Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter__tp_iternext,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)NULL,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};




static void
Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter__tp_clear(Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

}


static int
Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter__tp_traverse(Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter *self, visitproc visit, void *arg)
{
    Py_VISIT((PyObject *) self->container);
    return 0;
}


static void
_wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt____tp_dealloc(Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__ *self)
{
    delete self->obj;
    self->obj = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static void
_wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter__tp_dealloc(Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter *self)
{
    Py_CLEAR(self->container);
    delete self->iterator;
    self->iterator = NULL;

    self->ob_type->tp_free((PyObject*)self);
}


static PyObject*
_wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt____tp_iter(Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__ *self)
{
    Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter *iter = PyObject_GC_New(Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter, &Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter_Type);
    Py_INCREF(self);
    iter->container = self;
    iter->iterator = new std::vector<std::vector<std::vector<int> > >::iterator(self->obj->begin());
    return (PyObject*) iter;
}


static PyObject*
_wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter__tp_iter(Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter *self)
{
    Py_INCREF(self);
    return (PyObject*) self;
}

static PyObject* _wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter__tp_iternext(Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter *self)
{
    PyObject *py_retval;
    std::vector<std::vector<std::vector<int> > >::iterator iter;
    Pystd__vector__lt__std__vector__lt__int__gt_____gt__ *py_std__vector__lt__std__vector__lt__int__gt_____gt__;
    
    iter = *self->iterator;
    if (iter == self->container->obj->end()) {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    ++(*self->iterator);
    py_std__vector__lt__std__vector__lt__int__gt_____gt__ = PyObject_New(Pystd__vector__lt__std__vector__lt__int__gt_____gt__, &Pystd__vector__lt__std__vector__lt__int__gt_____gt___Type);
    py_std__vector__lt__std__vector__lt__int__gt_____gt__->obj = new std::vector<std::vector<int> >((*iter));
    py_retval = Py_BuildValue((char *) "N", py_std__vector__lt__std__vector__lt__int__gt_____gt__);
    return py_retval;
}

int _wrap_convert_py2c__std__vector__lt___std__vector__lt___std__vector__lt___int___gt_____gt_____gt__(PyObject *arg, std::vector<std::vector<std::vector<int> > > *container)
{
    if (PyObject_IsInstance(arg, (PyObject*) &Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt___Type)) {
        *container = *((Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__*)arg)->obj;
    } else if (PyList_Check(arg)) {
        container->clear();
        Py_ssize_t size = PyList_Size(arg);
        for (Py_ssize_t i = 0; i < size; i++) {
            std::vector< std::vector< int > > item;
            if (!_wrap_convert_py2c__std__vector__lt___std__vector__lt___int___gt_____gt__(PyList_GET_ITEM(arg, i), &item)) {
                return 0;
            }
            container->push_back(item);
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "parameter must be None, a vector_of_vector_of_vector_of_ints instance, or a list of std::vector< std::vector< int > >");
        return 0;
    }
    return 1;
}


static int
_wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt____tp_init(Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__ *self, PyObject *args, PyObject *kwargs)
{
    const char *keywords[] = {"arg", NULL};
    PyObject *arg = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, (char *) "|O", (char **) keywords, &arg)) {
        return -1;
    }

    self->obj = new std::vector<std::vector<std::vector<int> > >;

    if (arg == NULL)
        return 0;

    if (!_wrap_convert_py2c__std__vector__lt___std__vector__lt___std__vector__lt___int___gt_____gt_____gt__(arg, self->obj)) {
        delete self->obj;
        self->obj = NULL;
        return -1;
    }
    return 0;
}

PyTypeObject Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt___Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_vector_of_vector_of_ints",            /* tp_name */
    sizeof(Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt____tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)NULL,     /* tp_traverse */
    (inquiry)NULL,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt____tp_iter,          /* tp_iter */
    (iternextfunc)NULL,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)_wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt____tp_init,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};

PyTypeObject Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size */
    (char *) "PolytopeModules.vector_of_vector_of_vector_of_intsIter",            /* tp_name */
    sizeof(Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter),                  /* tp_basicsize */
    0,                                 /* tp_itemsize */
    /* methods */
    (destructor)_wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter__tp_dealloc,        /* tp_dealloc */
    (printfunc)0,                      /* tp_print */
    (getattrfunc)NULL,       /* tp_getattr */
    (setattrfunc)NULL,       /* tp_setattr */
    (cmpfunc)NULL,           /* tp_compare */
    (reprfunc)NULL,             /* tp_repr */
    (PyNumberMethods*)NULL,     /* tp_as_number */
    (PySequenceMethods*)NULL, /* tp_as_sequence */
    (PyMappingMethods*)NULL,   /* tp_as_mapping */
    (hashfunc)NULL,             /* tp_hash */
    (ternaryfunc)NULL,          /* tp_call */
    (reprfunc)NULL,              /* tp_str */
    (getattrofunc)NULL,     /* tp_getattro */
    (setattrofunc)NULL,     /* tp_setattro */
    (PyBufferProcs*)NULL,  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT|Py_TPFLAGS_HAVE_GC,                      /* tp_flags */
    NULL,                        /* Documentation string */
    (traverseproc)Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter__tp_traverse,     /* tp_traverse */
    (inquiry)Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter__tp_clear,             /* tp_clear */
    (richcmpfunc)NULL,   /* tp_richcompare */
    0,             /* tp_weaklistoffset */
    (getiterfunc)_wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter__tp_iter,          /* tp_iter */
    (iternextfunc)_wrap_Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter__tp_iternext,     /* tp_iternext */
    (struct PyMethodDef*)NULL, /* tp_methods */
    (struct PyMemberDef*)0,              /* tp_members */
    NULL,                     /* tp_getset */
    NULL,                              /* tp_base */
    NULL,                              /* tp_dict */
    (descrgetfunc)NULL,    /* tp_descr_get */
    (descrsetfunc)NULL,    /* tp_descr_set */
    0,                 /* tp_dictoffset */
    (initproc)NULL,             /* tp_init */
    (allocfunc)PyType_GenericAlloc,           /* tp_alloc */
    (newfunc)PyType_GenericNew,               /* tp_new */
    (freefunc)0,             /* tp_free */
    (inquiry)NULL,             /* tp_is_gc */
    NULL,                              /* tp_bases */
    NULL,                              /* tp_mro */
    NULL,                              /* tp_cache */
    NULL,                              /* tp_subclasses */
    NULL,                              /* tp_weaklist */
    (destructor) NULL                  /* tp_del */
};




PyMODINIT_FUNC
#if defined(__GNUC__) && __GNUC__ >= 4
__attribute__ ((visibility("default")))
#endif
initPolytopeModules(void)
{
    PyObject *m;
    PyObject *submodule;
    m = Py_InitModule3((char *) "PolytopeModules", PolytopeModules_functions, NULL);
    if (m == NULL) {
        return;
    }
    /* Register the 'std::set<unsigned>' class */
    if (PyType_Ready(&Pystd__set__lt__unsigned__gt___Type)) {
        return;
    }
    if (PyType_Ready(&Pystd__set__lt__unsigned__gt__Iter_Type)) {
        return;
    }
    PyModule_AddObject(m, (char *) "set_of_uints", (PyObject *) &Pystd__set__lt__unsigned__gt___Type);
    PyModule_AddObject(m, (char *) "set_of_uintsIter", (PyObject *) &Pystd__set__lt__unsigned__gt__Iter_Type);
    /* Register the 'std::vector<int>' class */
    if (PyType_Ready(&Pystd__vector__lt__int__gt___Type)) {
        return;
    }
    if (PyType_Ready(&Pystd__vector__lt__int__gt__Iter_Type)) {
        return;
    }
    PyModule_AddObject(m, (char *) "vector_of_ints", (PyObject *) &Pystd__vector__lt__int__gt___Type);
    PyModule_AddObject(m, (char *) "vector_of_intsIter", (PyObject *) &Pystd__vector__lt__int__gt__Iter_Type);
    /* Register the 'std::vector<unsigned>' class */
    if (PyType_Ready(&Pystd__vector__lt__unsigned__gt___Type)) {
        return;
    }
    if (PyType_Ready(&Pystd__vector__lt__unsigned__gt__Iter_Type)) {
        return;
    }
    PyModule_AddObject(m, (char *) "vector_of_uints", (PyObject *) &Pystd__vector__lt__unsigned__gt___Type);
    PyModule_AddObject(m, (char *) "vector_of_uintsIter", (PyObject *) &Pystd__vector__lt__unsigned__gt__Iter_Type);
    /* Register the 'std::vector<double>' class */
    if (PyType_Ready(&Pystd__vector__lt__double__gt___Type)) {
        return;
    }
    if (PyType_Ready(&Pystd__vector__lt__double__gt__Iter_Type)) {
        return;
    }
    PyModule_AddObject(m, (char *) "vector_of_doubles", (PyObject *) &Pystd__vector__lt__double__gt___Type);
    PyModule_AddObject(m, (char *) "vector_of_doublesIter", (PyObject *) &Pystd__vector__lt__double__gt__Iter_Type);
    /* Register the 'std::vector<std::vector<int> >' class */
    if (PyType_Ready(&Pystd__vector__lt__std__vector__lt__int__gt_____gt___Type)) {
        return;
    }
    if (PyType_Ready(&Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter_Type)) {
        return;
    }
    PyModule_AddObject(m, (char *) "vector_of_vector_of_ints", (PyObject *) &Pystd__vector__lt__std__vector__lt__int__gt_____gt___Type);
    PyModule_AddObject(m, (char *) "vector_of_vector_of_intsIter", (PyObject *) &Pystd__vector__lt__std__vector__lt__int__gt_____gt__Iter_Type);
    /* Register the 'std::vector<std::vector<unsigned> >' class */
    if (PyType_Ready(&Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt___Type)) {
        return;
    }
    if (PyType_Ready(&Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter_Type)) {
        return;
    }
    PyModule_AddObject(m, (char *) "vector_of_vector_of_uints", (PyObject *) &Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt___Type);
    PyModule_AddObject(m, (char *) "vector_of_vector_of_uintsIter", (PyObject *) &Pystd__vector__lt__std__vector__lt__unsigned__gt_____gt__Iter_Type);
    /* Register the 'std::vector<std::vector<double> >' class */
    if (PyType_Ready(&Pystd__vector__lt__std__vector__lt__double__gt_____gt___Type)) {
        return;
    }
    if (PyType_Ready(&Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter_Type)) {
        return;
    }
    PyModule_AddObject(m, (char *) "vector_of_vector_of_doubles", (PyObject *) &Pystd__vector__lt__std__vector__lt__double__gt_____gt___Type);
    PyModule_AddObject(m, (char *) "vector_of_vector_of_doublesIter", (PyObject *) &Pystd__vector__lt__std__vector__lt__double__gt_____gt__Iter_Type);
    /* Register the 'std::vector<std::set<unsigned> >' class */
    if (PyType_Ready(&Pystd__vector__lt__std__set__lt__unsigned__gt_____gt___Type)) {
        return;
    }
    if (PyType_Ready(&Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter_Type)) {
        return;
    }
    PyModule_AddObject(m, (char *) "vector_of_set_of_unsigned", (PyObject *) &Pystd__vector__lt__std__set__lt__unsigned__gt_____gt___Type);
    PyModule_AddObject(m, (char *) "vector_of_set_of_unsignedIter", (PyObject *) &Pystd__vector__lt__std__set__lt__unsigned__gt_____gt__Iter_Type);
    /* Register the 'std::vector<std::vector<std::vector<int> > >' class */
    if (PyType_Ready(&Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt___Type)) {
        return;
    }
    if (PyType_Ready(&Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter_Type)) {
        return;
    }
    PyModule_AddObject(m, (char *) "vector_of_vector_of_vector_of_ints", (PyObject *) &Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt___Type);
    PyModule_AddObject(m, (char *) "vector_of_vector_of_vector_of_intsIter", (PyObject *) &Pystd__vector__lt__std__vector__lt__std__vector__lt__int__gt_____gt_____gt__Iter_Type);
    submodule = initPolytopeModules_polytope();
    if (submodule == NULL) {
        return;
    }
    Py_INCREF(submodule);
    PyModule_AddObject(m, (char *) "polytope", submodule);
}
