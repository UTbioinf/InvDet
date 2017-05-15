from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t
import picos as pic
import networkx as nx
import cvxopt as cvx
import cvxopt.lapack
import numpy as np

cdef extern from "maxcut.h" namespace "loon":
    cdef cppclass CppMaxCut "loon::MaxCut":
        MaxCut() except+
        MaxCut(int threshold) except+
        void clear()
        void set_small_threshold(int threshold)
        void add_edge(int u, int v, int w)
        Py_ssize_t number_of_nodes() const
        Py_ssize_t number_of_edges() const
        int get_edge_u(Py_ssize_t i) const
        int get_edge_v(Py_ssize_t i) const
        int get_edge_w(Py_ssize_t i) const
        int get_node_rawid(Py_ssize_t i) const
        const vector[bool_t]& get_solution() const
        double get_value() const
        bool_t solve()
        bool_t is_bipartite()
        bool_t exact_algorithm()
        void max_spanning_tree()

cdef class MaxCut:
    cdef CppMaxCut _maxcut
    cdef float _value
    cdef list _solution
    cdef int _maxnode

    def __cinit__(self, threshold = 15, maxnode = 60):
        self._maxcut.set_small_threshold(threshold)
        self._value = float("-inf")
        self._maxnode = maxnode

    def clear(self):
        self._maxcut.clear()

    def set_small_threshold(self, int threshold):
        self._maxcut.set_small_threshold(threshold)

    def set_sdp_maxnode(self, int maxnode):
        self._maxnode = maxnode

    def add_edge(self, int u, int v, int w):
        self._maxcut.add_edge(u, v, w)

    def number_of_nodes(self):
        return self._maxcut.number_of_nodes()

    def number_of_edges(self):
        return self._maxcut.number_of_edges()

    def edge_iter(self):
        cdef int i
        cdef int n = self._maxcut.number_of_edges()
        cdef int u, v, w
        for i in xrange(n):
            u = self._maxcut.get_edge_u(i)
            v = self._maxcut.get_edge_v(i)
            w = self._maxcut.get_edge_w(i)
            yield (u, v, w)

    def edge_name_iter(self):
        cdef int i
        cdef int n = self._maxcut.number_of_edges()
        cdef int u, v, w
        for i in xrange(n):
            u = self._maxcut.get_edge_u(i)
            v = self._maxcut.get_edge_v(i)
            w = self._maxcut.get_edge_w(i)
            yield (self._maxcut.get_node_rawid(u), self._maxcut.get_node_rawid(v), w)
    
    def edges(self):
        cdef int i
        cdef int n = self._maxcut.number_of_edges()
        cdef int u, v, w
        cdef list ret = []
        for i in xrange(n):
            u = self._maxcut.get_edge_u(i)
            v = self._maxcut.get_edge_v(i)
            w = self._maxcut.get_edge_w(i)
            ret.append( (u, v, w) )
        return ret

    def edge_names(self):
        cdef int i
        cdef int n = self._maxcut.number_of_edges()
        cdef int u, v, w
        cdef list ret = []
        for i in xrange(n):
            u = self._maxcut.get_edge_u(i)
            v = self._maxcut.get_edge_v(i)
            w = self._maxcut.get_edge_w(i)
            ret.append( (self._maxcut.get_node_rawid(u), self._maxcut.get_node_rawid(v), w) )
        return ret

    def node_name_iter(self):
        cdef int i
        cdef int n = self._maxcut.number_of_nodes()
        for i in xrange(n):
            yield self._maxcut.get_node_rawid(i)

    def node_names(self):
        cdef int i
        cdef int n = self._maxcut.number_of_nodes()
        return [self._maxcut.get_node_rawid(i) for i in xrange(n)]

    def get_node_name(self, int i):
        return self._maxcut.get_node_rawid(i)

    def get_solution(self):
        return self._solution

    def get_value(self):
        return self._value

    def solve(self, int min_iteration = 100, int max_iteration = 10000, float min_ratio = 0.878, float max_ratio = 0.995):
        if self._maxcut.solve():
            self._solution = self._maxcut.get_solution()
        else:
            if self._maxnode >= self._maxcut.number_of_nodes():
                self.approx_878(min_iteration, max_iteration, min_ratio, max_ratio)
                if self._value < self._maxcut.get_value():
                    self._value = self._maxcut.get_value()
                    self._solution = self._maxcut.get_solution()

    def is_bipartite(self):
        if self._maxcut.is_bipartite():
            self._solution = self._maxcut.get_solution()
            self._value = float("-inf")
            return True
        return False

    def exact_algorithm(self):
        if self._maxcut.exact_algorithm():
            self._value = self._maxcut.get_value()
            self._solution = self._maxcut.get_solution()
            return True
        return False

    def max_spanning_tree(self):
        self._maxcut.max_spanning_tree()
        self._value = self._maxcut.get_value()
        self._solution = self._maxcut.get_solution()

    cpdef build_graph(self):
        G = nx.Graph()
        cdef int i
        cdef int n
        n = self._maxcut.number_of_edges()
        for i in xrange(n):
            G.add_edge(self._maxcut.get_edge_u(i),
                    self._maxcut.get_edge_v(i),
                    weight=self._maxcut.get_edge_w(i))
        return G

    cpdef approx_878(self, int min_iter = 100, int max_iter = 10000, float min_ratio = 0.878, float max_ratio = 0.995):
        cdef int N, i, j, cnt
        cdef float obj_sdp, obj, o
        G = self.build_graph()
        N = self.number_of_nodes()
        maxcut = pic.Problem()
        X = maxcut.add_variable('X', (N, N), 'symmetric')
        t_mat = 1/4.*nx.laplacian_matrix(G, nodelist=range(N))
        t_mat = t_mat.tocoo()
        L = pic.new_param('L', cvx.spmatrix(t_mat.data.tolist(), t_mat.row.tolist(), 
                t_mat.col.tolist(), size=t_mat.shape))
        maxcut.add_constraint(pic.tools.diag_vect(X) == 1)
        maxcut.add_constraint(X>>0)
        maxcut.set_objective('max', L|X)
        maxcut.solve(verbose = 0)
        V = X.value # Cholesky factorization
        cvxopt.lapack.potrf(V)
        for i in xrange(N):
            for j in xrange(i+1, N):
                V[i, j] = 0
        cnt = 0
        obj_sdp = maxcut.obj_value()
        obj = 0.0
        while (cnt < min_iter or (obj < min_ratio*obj_sdp and cnt < max_iter) ):
            r = cvx.normal(N, 1)
            x = cvx.matrix(np.sign(V * r))
            o = (x.T*L*x).value[0]
            if o > obj:
                x_cut = x
                obj = o
            cnt += 1
            if obj >= max_ratio * obj_sdp:
                break
        #print "Iterations: ", cnt
        self._solution = [(False if each == -1 else True) for each in x_cut]
        self._value = obj
        
