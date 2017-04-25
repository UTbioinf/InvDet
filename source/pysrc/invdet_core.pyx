from libcpp.string cimport string

cdef extern from "invdet_core.h" namespace "loon":
    cdef cppclass CppInvDector:
        void read(const string& fname)

        void gen_graphs(const string& fname)
        void gen_graphs(const string& fname, int min_cvg)
        void gen_graphs(const string& fname, int min_cvg, double min_cvg_percent)

cdef class InvDector:
    cdef CppInvDector _invdet

    cpdef read(self, str fname):
        self._invdet.read(<string>fname)

    cpdef gen_graphs(str fname, int min_cvg = 0, double min_cvg_percent = 0.0):
        self._invdet.gen_graphs(<string>fname, min_cvg, min_cvg_percent)

