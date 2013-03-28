#include "Python.h"
#include "math.h"
#include <stdio.h>
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"

#include <vista/VImage.h>

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

static PyObject *python2vista(PyObject *self, PyObject *arg)
{
    PyArrayObject *array;

    if (!PyArg_ParseTuple(arg, "O!", &PyArray_Type, &array))
        return NULL;
    if (array->nd != 2)
        return NULL;
    npy_intp *dim = array->dimensions;
    npy_intp *strd = array->strides;
    //unsigned int data[dim[0]][dim[2]];
    npy_intp strd0 = array->strides[0];
    npy_intp strd1 = array->strides[1];
    
    VImage image = VCreateImage (1, dim[0], dim[1], VLongRepn);
    printf("%i %i", dim[0],dim[1]);

    return PyString_FromString("sdkgfg");
}

static PyMethodDef VistaMethods[] = {
    {"python2vista", python2vista, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL} /*Sentinel*/
};

PyMODINIT_FUNC initvista()
{
/*
    PyObject *m;

    m = Py_InitModule("vista", VistaMethods);
    import_array();
    import_umath();
    return m;*/
    Py_InitModule3("vista", VistaMethods, "");
    import_array();
    import_umath();
}
