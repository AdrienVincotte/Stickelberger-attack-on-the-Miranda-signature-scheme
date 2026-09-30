# --- Paramètres du système ---
q = 2
m = 7
t = 1
la = 1
ls = 6

n = m
k = n - 2*t
assert m > k + la

load("config.sage")
load("utils.sage")
load("crypto.sage")
load("stickelberger.sage")

""" Key generation """

pk, sk = key_gen_chiffrement()
gamma = sk[1]
Gpub = pk

""" Truncation and reduction to a MinRank instance """

A = Matrix(Fq, kmasq, m*(k + la + 1))
V_mat = MatrixSpace(Fq, m, k + la + 1).random_element()
while V_mat.rank() != k + la + 1:
    V_mat = MatrixSpace(Fq, m, k + la + 1).random_element()

for i in range(kmasq):
    mati = pliage(Gpub[i], m, m)
    A[i] = (mati * V_mat).list()

ker = (A.transpose()).kernel()
base = ker.basis()
Mx = Matrix(Fq, m, k + la + 1)
for i in range(rho):
    Mx += x[i] * pliage(base[i], m, k + la + 1)
    
""" Construction of equations derived from minors """

equations = []
for i in range(m):
    for j in range(i + 1, m):
        for r_idx in range(k + la + 1):
            for s_idx in range(r_idx + 1, k + la + 1):
                equations.append(Mx[i, r_idx]*Mx[j, s_idx] - Mx[i, s_idx]*Mx[j, r_idx])

equations_Fq = [R(eq) for eq in equations]
equations_Fq.append(x[0] - 1)

I = R.ideal(equations_Fq)

solutions = stickelberger(I, Fq, Fqm, rho, max_essais=50)

bases = []
for i in range(len(solutions)):
    E_mat = matrix(Fq, m, k + la + 1)
    for j in range(rho):
        E_mat += solutions[i][j] * pliage(base[j], m, k + la + 1)
    bases.append(E_mat.transpose()[0])
    
#print(bases)

""" We verify that the basis gamma is indeed among those calculated. """
print(f"Successful attack: {gamma/gamma[0] in bases}")
