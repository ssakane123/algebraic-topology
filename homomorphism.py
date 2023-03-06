import numpy as np

class Homomorphism:
    def __init__(self, A, base_ring = "Z"):
        self.A = A.copy()
        self.L = np.eye(A.shape[0], dtype = int)
        self.R = np.eye(A.shape[1], dtype = int)
        self.base_ring = base_ring
    
    def __str__(self):
        return str(self.A)
    
    def add(self, i, j, c, row, column):
        if row:
            self.L[i] += self.L[j] * c
            self.A[i, :] += self.A[j, :] * c
        if column:
            self.R[:, i] += self.R[:, j] * c
            self.A[:, i] += self.A[:, j] * c
    
    def mul(self, i, c, row, column):
        if row:
            self.L[i] *= c
            self.A[i, :] *= c
        if column:
            self.R[:, i] *= c
            self.A[:, i] *= c
    
    def swap(self, i, j, row, column):
        if row:
            self.L[[i, j]] = self.L[[j, i]]
            self.A[[i, j]] = self.A[[j, i]]
        if column:
            self.R[:, [i, j]] = self.R[:, [j, i]]
            self.A[:, [i, j]] = self.A[:, [j, i]]
    
    def argmin(self, s):
        B = np.abs(self.A[s:, s:])
        C = np.nonzero(B)
        if C[0].size == 0:
            return (None, None)
        i, j = np.where(B == np.min(B[C]))
        return i[0] + s, j[0] + s
    
    def next_one(self, s):
        for i in range(s + 1, self.A.shape[0]):
            for j in range(s + 1, self.A.shape[1]):
                if self.A[i, j] % self.A[s, s]:
                    return (i, j)
        return (None, None)
    
    def is_nf(self, s):
        for i in range(s + 1, self.A.shape[0]):
            if self.A[i, s] != 0:
                return False
        for j in range(s + 1, self.A.shape[1]):
            if self.A[s, j] != 0:
                return False
        return True

    def _nf(self, s):
        if s == min(self.A.shape):
            return self.A
        
        i, j = self.argmin(s)
        if (i, j) == (None, None):
            return self.A
        
        self.swap(s, i, True, False)
        self.swap(s, j, False, True)

        for i in range(s + 1, self.A.shape[0]):
            if self.A[i, s] == 0:
                continue
            k = self.A[i, s] // self.A[s, s]
            self.add(i, s, -k, True, False)
        
        for j in range(s + 1, self.A.shape[1]):
            if self.A[s, j] == 0:
                continue
            k = self.A[s, j] // self.A[s, s]
            self.add(j, s, -k, False, True)
        
        if self.is_nf(s):
            i, j = self.next_one(s)
            if (i, j) != (None, None):
                self.add(s, i, 1, True, False)
                return self._nf(s)
            elif self.A[s, s] < 0:
                self.mul(s, -1, True, False)
            return self._nf(s + 1)
        else:
            return self._nf(s)
    
    def nf(self):
        A = self.A.copy()
        B = self._nf(0)
        self.A = A
        self.L = np.eye(A.shape[0], dtype = int)
        self.R = np.eye(A.shape[1], dtype = int)
        return B
