########### LR-Based IAG Decoder with Gaussian for weight r = \lfloor N(n-k)/(N+1) \rfloor ########

def random_small_vector_genenration(Extension, Length, Weight):
    B = matrix(Fqm.base_ring(), Weight, Extension, 0)
    while B.rank() != Weight:
        B = random_matrix(Fqm.base_ring(),Weight, Extension)
    C = matrix(Fqm.base_ring(), Length, Extension,0)
    while C.rank() != Weight:
        C = random_matrix(Fqm.base_ring(), Length, Weight) * B
    return vector(Fqm,[C[i] for i in range(Length)])

def random_support_genenration(t):
    B = matrix(Fqm.base_ring(),t,m,0)
    while B.rank() != t:
        B = matrix(Fqm.base_ring(),[vector(Fqm.random_element()) for i in range(t)])
    return vector(Fqm,[B[i] for i in range(t)])

def random_small_vector_genenration_with_support(n,support):
    t = len(list(support))
    C = matrix(Fqm.base_ring(),n,t,0)
    while C.rank() != t:
        C = matrix(Fqm.base_ring(),n,t,[Fqm.base_ring().random_element() for _ in range(n*t)]) 
    c = C * support
    return c

def Frob_Map(element,degree):
    a = element
    if degree == 0:
        Identity = Frob.inverse() * Frob
        a = Identity(a)
    if degree > 0:
        for i in range(degree):
            a = Frob(a)
    if degree < 0:
        InverseFrob = Frob.inverse()
        for i in range(-degree):
            a = InverseFrob(a)
    return a

def Frob_vector(vectors,degrees):
    lengths = len(list(vectors))
    Frob_vector = vector(Fqm, lengths, [Frob_Map(vectors[i], degrees) for i in range(lengths)])
    return Frob_vector

def Moore_matrix(vectors, nrows):
    List = []; ncolumns = len(list(vectors))
    for i in range(nrows):
        List.append(Frob_vector(vectors,i))
    return matrix(Fqm, nrows, ncolumns, List)

def test(totalltests):
    succ = 0
    failure = 0
    for npair in range(totalltests):
        f = [S.random_element(degree=(-1, k-1)) for _ in range(N)] # Random messages
        support = random_support_genenration(r)
        e = [random_small_vector_genenration_with_support(n, support) for _ in range(N)] # Random error
        y = [vector([f[i](g[j]) + e[i][j] for j in range(n)]) for i in range(N)] # Received word
        Y = [Moore_matrix(yi, r+1).transpose() for yi in y] # Moore matrices
        A1 = block_matrix(Fqm, N, 1, Y) # Build A1
        A = block_matrix(Fqm, 1, 2, [A1, A2]) # Total matrix
        Solution = A.right_kernel_matrix()[0].list() # Solve
        V = S(Solution[:r+1]) # V(x)
        # Recover all Ni
        start = r + 1
        success = True
        for i in range(N):
            coeff = Solution[start:start+k+r]
            Ni = S(coeff)
            fi, rem = Ni.left_quo_rem(-V)
            if (rem != 0) or (fi != f[i]):
                success = False
                break
            start += k+r

        if success:
            succ += 1
        else:
            failure += 1
    print("success/totalltests: %d/%d; success rate: %f" % (succ, totalltests, succ/totalltests))
    print("failure/totalltests: %d/%d; failure rate: %f" % (failure, totalltests, failure/totalltests))

# extension degree = m; code_length = n; generator_weight = t; code_dimension = k; rank weight = r; Interleaved order = N
# r = \lfloor N(n-k)/(N+1) \rfloor
# \gamma_q = 4

(q, m, n, t, k, r, N) = (2, 5, 5, 4, 2, 2, 3)  
# DFR: Simulated DFR: 0.319; Theoretical DFR:  4*2^{-2} = 1



Fqm = GF(q**m)
Frob = Fqm.frobenius_endomorphism()
S = OrePolynomialRing(Fqm, Frob, 'x')
# S.<x> = Fqm['x', Frob]
#print(S)

g = random_small_vector_genenration(m,n,min(m,n,t))  # Generator of AG codes
Z = Moore_matrix(g, k+r).transpose()
A2 = block_diagonal_matrix([Z] * N)



%time test(100000)
