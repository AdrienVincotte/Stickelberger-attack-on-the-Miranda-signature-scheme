""" Application of masking to the unfolded matrix """

def masquage(Gdep):

    Gdep_ls = Gdep[0:m*k - ls]
    for i in range(la):
        ligne = Mdeplie.random_element()
        while ligne in span(Gdep):
            ligne = Mdeplie.random_element()
        Gdep_ls = Gdep_ls.stack(Matrix(ligne))
    
    P = Brouillage.random_element()
    while not P.is_invertible():
        P = Brouillage.random_element()
    return P * Gdep_ls

""" Key generation """

def key_gen_chiffrement():

    g = C.random_element()
    while psiv_g(g, base_canon).rank() != n:
        g = C.random_element()
        
    Gabi = codes.GabidulinCode(Fqm, n, k, Fq, evaluation_points=g)
    G = Gabi.generator_matrix()
    gamma = tirage_base()
    G_deplie = vec_to_mat(G, gamma)
    
    sk = (G, gamma)
    pk = masquage(G_deplie)
    return (pk, sk)
