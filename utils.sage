""" Random sampling of an Fq-basis of Fqm """

def tirage_base():
    rep = MM.random_element()
    if rep.is_invertible():
        return BFqm([from_V(ligne) for ligne in rep])
    return tirage_base()

""" Folding of a vector of length mn into a matrix of dimension m*n """

def pliage(v, lignes, colonnes):
    assert len(v) == lignes * colonnes
    mat = matrix(Fq, lignes, colonnes)
    for i in range(lignes):
        mat[i] = v[i*colonnes:(i + 1)*colonnes]
    return mat

""" Transition matrix from the standard basis to a gamma basis (for horizontal vector) """

def matrice_passage(gamma):
    assert len(gamma) == m
    matp = [to_V(gamma[i]) for i in range(m)]
    return MM(matp).inverse()

""" Coordinates of an element of Fqm in an Fq-basis gamma """

def psix_g(x_elem, gamma):
    return to_V(x_elem) * matrice_passage(gamma)

"""Extends a vector in Fqm of length n into an m × n matrix"""

def psiv_g(v, gamma):
    assert n == len(v)
    psig = [psix_g(v[i], gamma) for i in range(n)]
    return Mtransp(psig).transpose()

""" Matrix of multiplication by 'a' in the vertical canonical basis """

def multiplication_fqm():
    Pmin = a.minimal_polynomial()
    Mult_a = matrix(Fq, m, m)
    for i in range(m - 1):
        Mult_a[i + 1, i] = 1
    for j in range(m):
        Mult_a[j, m - 1] = Pmin[j]
    return Mult_a

""" Transforms an Fqm-linear vector code into an extended matrix code in the gamma basis """

def vec_to_mat(G, gamma):
    Mult_a = multiplication_fqm()
    Gdeplie = matrix(Fq, k*m, n*m)
    P = matrice_passage(gamma).transpose()
    Mult_a_gamma = P * Mult_a * (P^(-1))
    i = 0
    for ligne in G:
        L = psiv_g(ligne, gamma)
        for j in range(m):
            Gdeplie[i] = [val for row in L.rows() for val in row]
            L = Mult_a_gamma * L
            i += 1
    P = Melange.random_element()
    while not P.is_invertible():
        P = Melange.random_element()
    return P * Gdeplie
