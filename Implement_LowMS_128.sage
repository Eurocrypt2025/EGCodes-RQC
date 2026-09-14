Implement LowMS-128

def random_small_vector_genenration(Extension, Length, Weight):
    B = matrix(Fqm.base_ring(), Weight, Extension, 0)
    while B.rank() != Weight:
        B = random_matrix(Fqm.base_ring(),Weight, Extension)
    C = matrix(Fqm.base_ring(), Length, Extension,0)
    while C.rank() != Weight:
        C = random_matrix(Fqm.base_ring(), Length, Weight) * B
    return vector(Fqm,[C[i] for i in range(Length)])

def random_small_space_gen(t):
    B = matrix(Fqm.base_ring(),t,m,0)
    while B.rank() != t:
        B = matrix(Fqm.base_ring(),[vector(Fqm.random_element()) for i in range(t)])
    return B.row_space()

def homogenous_matrix_gen(k,n,d): # generate k * n homogenous matrix of dimension d
    F = random_small_space_gen(d)
    H_list = [Fqm(list(F.random_element())) for _ in range(n*k)]
    return matrix(Fqm,k,n,H_list)


def matrix_to_basis_list(B):
    b = vector(B)
    c = B.nrows() * B.ncols()
    BB =  matrix(Fqm.base_ring(),c,m,[vector(b[j]) for j in range(c)]).row_space().basis_matrix()
    return [Fqm(x) for x in BB.rows()]

def vector_matrix(small_vector): # 向量、列表、多项式、转化成矩阵
    length = len(list(small_vector))
    return matrix(Fqm.base_ring(),length,m,[vector(small_vector[j]) for j in range(length)])

def vector_space(small_vector): # 求向量、列表、多项式决定的行空间
    return vector_matrix(small_vector).row_space()

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


# Key Generation
def LowMS_KGen(q,m,n,k,r,w):
    S = random_matrix(Fqm, n-k, n-k)
    while S.inverse() == 0:
        S = random_matrix(Fqm, n-k, n-k)
    g = random_small_vector_genenration(m,n,min(m,n))  # Generator of AG codes
    Gabi = matrix(Fqm,k,n,[g[i]**(q**j) for j in range(k) for i in range(n)]) # 生成矩阵
    H = Gabi.right_kernel_matrix()  # 校验矩阵
    P = homogenous_matrix_gen(n, n, w)
    while P.inverse() == 0:
        P = homogenous_matrix_gen(n, n, w)
    HH = S*H*P
    return S.inverse(), H, P.inverse(), g, HH



# Encapsulation
def LowMS_Encap(Public_Key):  
    HH = Public_Key
    e = random_small_vector_genenration(m, 6*n, r)   
    e1 = e[:n]; e2 = e[n:2*n]; e3 = e[2*n:3*n]; e4 = e[3*n:4*n]; e5 = e[4*n:5*n]; e6 = e[5*n:6*n]
    s1 = HH*e1; s2 = HH*e2; s3 = HH*e3; s4 = HH*e4; s5 = HH*e5; s6 = HH*e6
    C = matrix(N, n-k, [s1, s2, s3, s4, s5, s6]).transpose()
    EE = matrix(N, n, [e1, e2, e3, e4, e5, e6]).transpose()
    E = vector_space(matrix_to_basis_list(EE))
    E_string = "".join(map(str,list(E.basis_matrix()))) # list(E.basis_matrix()) can be replace by E.
    return C, hashlib.sha256(E_string.encode()).hexdigest()


# Decapsulation 
def LowMS_Decap(Private_Key, Ciphertext):  
    S = Private_Key[0]; H = Private_Key[1]; P = Private_Key[2]; g = Private_Key[3]
    C = Ciphertext
    SC = (S * C).transpose()

    ##### decoding IAG codes by Gaussian elimination #####
    s1 = SC[0]; s2 = SC[1]; s3 = SC[2]; s4 = SC[3]; s5 = SC[4]; s6 = SC[5]
    y1 = H.solve_right(s1); y2 = H.solve_right(s2); y3 = H.solve_right(s3)
    y4 = H.solve_right(s4); y5 = H.solve_right(s5); y6 = H.solve_right(s6)
    r = 21 
    Z = Moore_matrix(g, k + r).transpose()
    A2 = block_diagonal_matrix(Z, Z, Z, Z, Z,  Z) 

    Y1 = Moore_matrix(y1, r+1).transpose(); Y2 = Moore_matrix(y2, r+1).transpose()
    Y3 = Moore_matrix(y3, r+1).transpose(); Y4 = Moore_matrix(y4, r+1).transpose()
    Y5 = Moore_matrix(y5, r+1).transpose(); Y6 = Moore_matrix(y6, r+1).transpose()
    A1 = block_matrix(N, 1, [Y1, Y2, Y3, Y4, Y5, Y6]) 
    A = block_matrix(Fqm, 1, 2, [A1, A2])

    Solution = A.right_kernel_matrix()[0].list()  
    V = OPR(Solution[0: r+1])
    N1_vector = vector(Solution[r+1: r+1 + k + r])
    N2_vector = vector(Solution[r+1 + k + r: r+1 + 2*(k + r)])
    N3_vector = vector(Solution[r+1 + 2*(k + r): r+1 + 3*(k + r)])
    N4_vector = vector(Solution[r+1 + 3*(k + r): r+1 + 4*(k + r)])
    N5_vector = vector(Solution[r+1 + 4*(k + r): r+1 + 5*(k + r)])
    N6_vector = vector(Solution[r+1 + 5*(k + r): r+1 + 6*(k + r)])

    N1 = OPR(list(N1_vector)); N2 = OPR(list(N2_vector)); N3 = OPR(list(N3_vector))
    N4 = OPR(list(N4_vector)); N5 = OPR(list(N5_vector)); N6 = OPR(list(N6_vector))
    ff1, re = N1.left_quo_rem(-V); ff2, re = N2.left_quo_rem(-V); ff3, re = N3.left_quo_rem(-V)
    ff4, re = N4.left_quo_rem(-V); ff5, re = N5.left_quo_rem(-V); ff6, re = N6.left_quo_rem(-V)

    ee1 = vector([y1[i] - ff1(g[i]) for i in range(n)]) 
    ee2 = vector([y2[i] - ff2(g[i]) for i in range(n)]) 
    ee3 = vector([y3[i] - ff3(g[i]) for i in range(n)]) 
    ee4 = vector([y4[i] - ff4(g[i]) for i in range(n)]) 
    ee5 = vector([y5[i] - ff5(g[i]) for i in range(n)]) 
    ee6 = vector([y6[i] - ff6(g[i]) for i in range(n)]) 

    PE = matrix(N, n,[ee1, ee2, ee3, ee4, ee5, ee6]).transpose()  # (PE = P*E)

    E_test = vector_space(matrix_to_basis_list(P * PE))
    E_string = "".join(map(str,list(E_test.basis_matrix()))) # list(E.basis_matrix()) can be replace by E.
    return hashlib.sha256(E_string.encode()).hexdigest()


(q, m, n, k, r, w, N) = (2, 61, 50, 25, 7, 3, 6)   # 128

import hashlib
Fqm = GF(q**m)
Frob = Fqm.frobenius_endomorphism()
OPR = OrePolynomialRing(Fqm, Frob, 'x')


%time Key_List = LowMS_KGen(q,m,n,k,r,w)

Public_Key = Key_List[4]
%time Encapsulation_Key_List = LowMS_Encap(Public_Key)

Private_Key = Key_List[0:4]; Ciphertext = Encapsulation_Key_List[0]
%time Key_test = LowMS_Decap(Private_Key, Ciphertext)

# Check Correctness
Key_test == Encapsulation_Key_List[1]
