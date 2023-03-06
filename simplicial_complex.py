import numpy as np
from collections import defaultdict
from homomorphism import Homomorphism

class SimplicialComplex:
    def __init__(self, K, base_ring = "Z"):
        self.__base_ring = base_ring
        self.__dims = defaultdict(list)
        self.__dim = 0
        for v in K:
            v = tuple(sorted(v))
            self.__dims[len(v) - 1].append(v)
            self.__dim = max(self.__dim, len(v) - 1)
        self.__chain = {q: len(self.__dims[q]) for q in range(-1, self.__dim + 1)}
        self.__boundary_homomorphism = {q: 0 for q in range(self.__dim + 1)}
        self.__cycle = {q: 0 for q in range(self.__dim + 1)}
        self.__boundary = {q: defaultdict(int) for q in range(self.__dim + 2)}
        self.__homology_group = {q: defaultdict(int) for q in range(self.__dim + 1)}
        self.__boundaryetti_number = {q: 0 for q in range(self.__dim + 1)}
        for q in range(self.__dim + 1):
            self._boundary_homomorphism(q)
        for q in range(self.__dim + 1):
            self._cycle_and_boundary(q)
        self.__boundary[self.__dim + 1] = {0: 0}
        for q in range(self.__dim + 1):
            self._homology_group(q)
            self._betti_number(q)
    
    def _boundary_homomorphism(self, q):
        if q <= 0 or q > self.__dim:
            self.__boundary_homomorphism[q] = 0
        else:
            m, n = self.__chain[q - 1], self.__chain[q]
            A = np.zeros((m, n), dtype = int)
            M = {v: j for j, v in enumerate(self.__dims[q])}
            N = {v: j for j, v in enumerate(self.__dims[q - 1])}
            for v in self.__dims[q]:
                for i in range(q + 1):
                    A[N[v[:i] + v[i + 1:]], M[v]] = 1 if i % 2 == 0 else -1
            D = Homomorphism(A).nf()
            self.__boundary_homomorphism[q] = np.array(D)

    def _cycle_and_boundary(self, q):
        if q <= -1 or q > self.__dim:
            self.__cycle[q] = 0
            self.__boundary[q] = {0: 0}
        elif q == 0:
            self.__cycle[q] = self.__chain[q]
            self.__boundary[q] = {0: 0}
        else:
            D = np.diag(self.__boundary_homomorphism[q])
            for i in np.nonzero(D)[0]:
                a = self.__boundary_homomorphism[q][i, i]
                if a == 1:
                    self.__boundary[q][a - 1] += 1
                else:
                    self.__boundary[q][a] += 1
            self.__cycle[q] = self.__chain[q] - sum(self.__boundary[q].values())

    def _homology_group(self, q):
        if q <= -1 or q > self.__dim:
            self.__homology_group[q] = 0
        else:
            self.__homology_group[q][0] = self.__cycle[q] - sum(self.__boundary[q + 1].values())
            for e, v in self.__boundary[q + 1].items():
                if e == 0:
                    continue
                self.__homology_group[q][e] = v

    def _betti_number(self, q):
        if self.__homology_group[q] == -1:
            self.get_homology_groups(q)
        self.__boundaryetti_number[q] = self.__homology_group[q][0]
    
    def _ring(self, q, e=0, quotient=False):
        digit = {0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}
        s = ""
        if q == 0:
            s = "0"
        elif q == 1:
            if e == 0:
                s = self.__base_ring
            else:
                s = "{}{}".format(e, self.__base_ring) if not quotient else "{}/{}".format(self.__base_ring, e)
        else:
            p = "".join(digit[int(t)] for t in map(int, str(q)))
            if e == 0:
                s = self.__base_ring + p
            else:
                s = "({}{})".format(e, self.__base_ring) + p if not quotient else "({}/{})".format(self.__base_ring, e) + p
        return s
    
    def chain(self, q=None):
        res = []
        if q is not None:
            res.append("C[{}] = {}".format(q, self._ring(self.__chain[q])))
        else:
            for q in range(self.__dim + 1):
                res.append("C[{}] = {}".format(q, self._ring(self.__chain[q])))
        print(*res, sep="\n")
    
    def boundary_homomorphism(self, q=None):
        res = []
        if q is not None:
            res[q] = "d[{}]: {} -> {}".format(q, self._ring(self.__chain[q]), self._ring(self.__chain[q - 1]))
            res[q] += "\n" + str(self.__boundary_homomorphism[q])
        else:
            for q in range(self.__dim + 1):
                res.append("d[{}]: {} -> {}".format(q, self._ring(self.__chain[q]), self._ring(self.__chain[q - 1])))
                res.append(str(self.__boundary_homomorphism[q]))
        print(*res, sep="\n")
    
    def cycle(self, q=None):
        res = []
        if q is not None:
            res.append("Z[{}]: {}".format(q, self._ring(self.__cycle[q])))
        else:
            for q in range(self.__dim + 1):
                res.append("Z[{}]: {}".format(q, self._ring(self.__cycle[q])))
        print(*res, sep="\n")
    
    def boundary(self, q=None):
        res = []
        if q is not None:
            tmp = []
            for e, v in self.__boundary[q].items():
                tmp.append(self._ring(v, e))
            res.append("B[{}]: ".format(q) + " ⊕ ".join(tmp))
        else:
            for q in range(self.__dim + 1):
                tmp = []
                for e, v in self.__boundary[q].items():
                    tmp.append(self._ring(v, e))
                res.append("B[{}]: ".format(q) + " ⊕ ".join(tmp))
        print(*res, sep="\n")

    def homology_group(self, q=None):
        res = []
        if q is not None:
            tmp = []
            for e, v in self.__homology_group[q].items():
                tmp.append(self._ring(v, e, True))
            res.append("H[{}]: ".format(q) + " ⊕ ".join(tmp))
        else:
            for q in range(self.__dim + 1):
                tmp = []
                for e, v in self.__homology_group[q].items():
                    tmp.append(self._ring(v, e, True))
                res.append("H[{}]: ".format(q) + " ⊕ ".join(tmp))
        print(*res, sep="\n")

    def betti_number(self, q=None):
        res = []
        if q is not None:
            res.append("b[{}]: {}".format(q, self.__boundaryetti_number[q]))
        else:
            for q in range(self.__dim + 1):
                res.append("b[{}]: {}".format(q, self.__boundaryetti_number[q]))
        print(*res, sep="\n")
    
    def __str__(self):
        return str(self.__dims)