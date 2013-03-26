#include "Python.h"
#include "math.h"
#include <stdio.h>
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include "numpy/halffloat.h"

/*
 * vista.c
 * This is the C code for creating your own
 * Numpy ufunc.
 *
 * In this code we only define the ufunc for
 * a single dtype. The computations that must
 * be replaced to create a ufunc for
 * a different function are marked with BEGIN
 * and END.
 *
 * 'Extending and Embedding' and 'Python/C API' at
 * docs.python.org .
 *
 */

static PyMethodDef VistaMethods[] = {
        {NULL, NULL, 0, NULL}
};

/* The loop definitions must precede the PyMODINIT_FUNC. */

static void python2vista(char **args, npy_intp *dimensions,
                  npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *in = args[0], *out = args[1];
    npy_intp in_step = steps[0], out_step = steps[1];

    /*BEGIN main part*/
    /*FILE *f;
    f = fopen('tt', 'w');
    fprintf(f, "tt");*/
    printf("args %f",*(double *)in);
    printf("dimension %i",*dimensions);
    /*END main part*/
}

/*This gives pointers to the above functions*/
PyUFuncGenericFunction funcs[1] = {&python2vista};

/* These are the input and return dtype.*/
static char types[2] = {NPY_DOUBLE, NPY_DOUBLE};

static void *data[1] = {NULL};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "vista",
    NULL,
    -1,
    VistaMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit_vista(void)
{
    PyObject *m, *python2vista, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    python2vista = PyUFunc_FromFuncAndData(funcs, data, types, 1, 1, 1,
                                           PyUFunc_None, "python2vista",
                                           "", 0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "python2vista", python2vista);
    Py_DECREF(python2vista);

    return m;
}
#else
PyMODINIT_FUNC initvista(void)
{
    PyObject *m, *python2vista, *d;

    m = Py_InitModule("vista", VistaMethods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    python2vista = PyUFunc_FromFuncAndData(funcs, data, types, 1, 1, 1,
                                           PyUFunc_None, "python2vista",
                                           "", 0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "python2vista", python2vista);
    Py_DECREF(python2vista);
}
#endif
