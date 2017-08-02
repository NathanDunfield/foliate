from libcpp.string cimport string
from libcpp.functional cimport function
import os, sys, re, tempfile

cdef extern from "engine.h" namespace "regina":
    cdef int versionMajor()
    cdef int versionMinor()

cdef extern from "maths/integer.h" namespace "regina":
    cdef cppclass LargeInteger:
        long longValue()

cdef extern from "triangulation/dim3.h" namespace "regina":
    cdef cppclass NTriangulation:
        NTriangulation()
        NTriangulation(string)
        unsigned long countTetrahedra()
        bint intelligentSimplify()
        bint simplifyToLocalMinimum(bint perform)
        bint retriangulate(int, unsigned int, void*, function[bint(NTriangulation&)])
        string isoSig()

cdef bint action(NTriangulation& triangulation):
    if not triangulation.simplifyToLocalMinimum(False):
        print(triangulation.isoSig())
        
def version():
    return "%d.%d" % (versionMajor(), versionMinor())

cdef class Triangulation:
    """
    A Regina triangulation of a 3-manifold, which can be specified
    either by name of a triangulation file in SnapPea format or by giving
    a SnapPy triangulation.
    """
    cdef NTriangulation* triangulation

    def __cinit__(self, spec=None):
        self.triangulation = new NTriangulation(spec)

    def __dealloc__(self):
        del self.triangulation
            
    def num_tetrahedra(self):
        return self.triangulation.countTetrahedra()

    def simplify(self):
        cdef bint did_simplify = self.ntriangulation.intelligentSimplify()
        return did_simplify

    def isosig(self):
        return self.triangulation.isoSig()

    def retriangulate(self, height, threads=1):
        cdef int nthreads
        nthreads = threads
        cdef function[bint(NTriangulation&)] *callback
        callback = new function[bint(NTriangulation&)](action)
        self.triangulation.retriangulate(height, nthreads, NULL, callback[0])
