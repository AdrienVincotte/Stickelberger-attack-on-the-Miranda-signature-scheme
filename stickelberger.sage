def stickelberger(I, Fq, Fqm, rho, max_essais=50):
    R = I.ring()
    x = R.gens()

    # 1. Normal basis of the quotient algebra A' = Fq[x]/I
    B = I.normal_basis()
    dim_A = len(B)
    t2 = time.time()

    # 2. Projection matrix C (rho x dim_A) on Fq :
    # Each variable x[idx_var] modulo I is expressed in the basis B.
    C = matrix(Fq, rho, dim_A)
    for idx_var in range(rho):
        poly_reduit = x[idx_var].reduce(I)
        for idx_b, monom in enumerate(B):
            C[idx_var, idx_b] = poly_reduit.monomial_coefficient(monom)

    C_ext = C.change_ring(Fqm)
    I_dim = identity_matrix(Fqm, dim_A)
    
    # 3. Random sampling of coefficients from Fq (x0 is excluded as it is fixed at 1)
    for essai in range(1, max_essais + 1):
        coeffs = [Fq(0)] + [Fq.random_element() for _ in range(1, rho)]
        if all(c == 0 for c in coeffs):
            continue

        L = sum(coeffs[idx_var] * x[idx_var] for idx_var in range(rho))

        # 4. Construction of matrix associated to L
        M_L = matrix(Fq, dim_A, dim_A)
        for j_b, b in enumerate(B):
            poly_reduit = (L * b).reduce(I)
            for i_b, monom in enumerate(B):
                M_L[i_b, j_b] = poly_reduit.monomial_coefficient(monom)

        # 5. Transition to Fqm and transposition to work in the dual space.
        M_L_ext = M_L.change_ring(Fqm)
        M_L_T = M_L_ext.transpose()

        # 6. Search for eigenvalues ​​in Fqm
        P_char = M_L.change_ring(Fqm).characteristic_polynomial()
        racines_L = P_char.roots()

        # Verification that L separates all the roots in Fqm
        if len(racines_L) != dim_A or any(mult != 1 for _, mult in racines_L):
            #print(f"[-] Essai {essai} : forme non separatrice ou racines non simples.")
            continue

        separatrice = True
        solutions_essai = []

        # 7. Simultaneous reconstruction of coordinates for each eigenspace
        for val_L, _ in racines_L:
            noyau = (M_L_T - val_L * I_dim).right_kernel()

            if noyau.dimension() != 1:
                separatrice = False
                break

            v = noyau.basis()[0]

            # B[0] corresponds to the monomial 1; v[0] must be non-zero for normalization.
            if v[0] == 0:
                separatrice = False
                break

            v_norm = v / v[0]

            # 8. Reconstruction of the complete solution vector: sol = C * v_norm
            sol = C_ext * v_norm
            sol = sol/sol[0]
            solutions_essai.append(sol)

        if separatrice and len(solutions_essai) == dim_A:
            #print(f"[+] Forme séparatrice valide trouvée à l'essai {essai}.")
            return solutions_essai

    raise RuntimeError("Unable to find a separating shape after several attempts.")
