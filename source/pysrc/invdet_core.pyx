from libcpp.string cimport string

cdef extern from "invdet_core.h" namespace "loon":
    cdef cppclass CppInvDector "loon::InvDector":
        void read(const string& fname) except +RuntimeError

        void gen_graphs(const string& fname, int min_cvg, double min_cvg_percent, int min_overlap) except +RuntimeError
        void gen_graphs(const string& fname, int min_cvg, double min_cvg_percent, int min_overlap, const string& nucmer_prefix) except +RuntimeError
        
        void report_inversions(const string& graph_fname, const string& maxcut_fname,
                const string& inversion_fname) except +RuntimeError

cdef class InvDector:
    cdef CppInvDector _invdet

    def read(self, str fname):
        self._invdet.read(<string>fname)

    def gen_graphs(self, str fname, int min_cvg = 0, double min_cvg_percent = 0.0, int min_overlap=0):
        self._invdet.gen_graphs(<string>fname, min_cvg, min_cvg_percent, min_overlap)

    def gen_graphs_ignore_inverted_repeats(self, str fname, str nucmer_prefix, int min_cvg = 0, double min_cvg_percent = 0.0, int min_overlap=0):
        self._invdet.gen_graphs(<string>fname, min_cvg, min_cvg_percent, min_overlap, <string>nucmer_prefix)

    def report_inversions(self, str graph_fname, str maxcut_fname, str inversion_fname):
        self._invdet.report_inversions( <string>graph_fname, <string>maxcut_fname, <string>inversion_fname )

