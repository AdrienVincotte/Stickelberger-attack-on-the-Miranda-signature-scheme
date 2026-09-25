def tirage_base():
    """Tirage aléatoire d'une Fq-base de Fqm."""
    rep = MM.random_element()
    if rep.is_invertible():
        return BFqm([from_V(ligne) for ligne in rep])
    return tirage_base()


def pliage(v, lignes, colonnes):
    """Pliage d'un vecteur de longueur mn en matrice de taille m*n."""
    assert len(v) == lignes * colonnes
    mat = matrix(Fq, lignes, colonnes)
    for i in range(lignes):
        mat[i] = v[i*colonnes:(i + 1)*colonnes]
    return mat


def matrice_passage(gamma):
    """Matrice de passage de la base canonique vers une base gamma (vecteur horizontal)."""
    assert len(gamma) == m
    matp = [to_V(gamma[i]) for i in range(m)]
    return MM(matp).inverse()


def psix_g(x_elem, gamma):
    """Écrit un élément de Fqm dans une Fq-base gamma."""
    return to_V(x_elem) * matrice_passage(gamma)


def psiv_g(v, gamma):
    """Transforme un vecteur en matrice dont les coordonnées sont en colonnes."""
    assert n == len(v)
    psig = [psix_g(v[i], gamma) for i in range(n)]
    return Mtransp(psig).transpose()


def multiplication_fqm():
    """Matrice de multiplication par 'a' dans la base canonique verticale."""
    Pmin = a.minimal_polynomial()
    Mult_a = matrix(Fq, m, m)
    for i in range(m - 1):
        Mult_a[i + 1, i] = 1
    for j in range(m):
        Mult_a[j, m - 1] = Pmin[j]
    return Mult_a


def vec_to_mat(G, gamma):
    """Transforme un code vectoriel Fqm-linéaire en code matriciel étendu en base gamma."""
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