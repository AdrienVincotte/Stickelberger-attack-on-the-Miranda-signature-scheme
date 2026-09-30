import time
import random

""" Base fields and extensions """

Fqm.<a> = GF(q^m)
Fq = Fqm.base_ring()
V, from_V, to_V = Fqm.vector_space(GF(q), map=True)
Pmin = a.minimal_polynomial()

""" Vector and matrix spaces """

C = VectorSpace(Fqm, n)
BFqm = VectorSpace(Fqm, m)
M = MatrixSpace(Fq, m, n)
Mtransp = MatrixSpace(Fq, n, m)
Mdeplie = VectorSpace(Fq, m*n)
MM = MatrixSpace(Fq, m, m)
kmasq = k*m - ls + la
Brouillage = MatrixSpace(Fq, kmasq, kmasq)

base_canon = BFqm([a^i for i in range(m)])
Melange = MatrixSpace(Fq, k*m, k*m)

""" Polynomial ring enabling to solve the MinRank instance """

rho = (m-1)*(la+1) + ls + 1
variables = [f"x{i}" for i in range(rho)]
R = PolynomialRing(Fq, names=variables)
R.inject_variables(verbose=False)
x = R.gens()
