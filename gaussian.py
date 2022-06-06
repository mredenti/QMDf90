import numpy as np

class Gaussian():

    def __init__(self, a, q, p, C, eps, s=0, state=None):

        self.q = q
        self.p = p
        self.C = C
        self.eps = eps
        self.q = q
        self.s = s 
        self.dim =1 
        self.state = state
        self.a = a

    def psi(self, x):

        detC = np.linalg.det(self.C.imag) if self.dim != 1 else self.C.imag

        return np.exp(1j/self.eps * self.s.real) * self.a * (np.pi * self.eps)**(-self.dim / 4) * detC**(0.25) * np.exp(1j / 2 / self.eps * (x - self.q).T * self.C * (x - self.q) + 1j / self.eps*self.p.T*(x-self.q))

    def l2norm(self):

        x, dx = np.linspace(self.q - 4, self.q + 4, 10**5, retstep=True)

        return np.sum(np.abs(self.psi(x))**2) * dx  