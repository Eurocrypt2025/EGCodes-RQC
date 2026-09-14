###############  Decoding AG Codes with Improved Gaussian Elimination algorithm ###############

def random_small_vector_genenration(Extension, Length, Weight):
    B = matrix(Fqm.base_ring(), Weight, Extension, 0)
    while B.rank() != Weight:
        B = random_matrix(Fqm.base_ring(),Weight, Extension)
    C = matrix(Fqm.base_ring(), Length, Extension,0)
    while C.rank() != Weight:
        C = random_matrix(Fqm.base_ring(), Length, Weight) * B
    return vector(Fqm,[C[i] for i in range(Length)])

def Decoding_AGCodes(Noisy_Word, r): 
    y_monomials = [Noisy_Word[i]**(q**j) for i in range(n) for j in range(r+1)] 
    T = matrix(Fqm, n, r+1, y_monomials) 
    T1 = T.matrix_from_rows(range(k+r))
    T2 = T.matrix_from_rows(range(k+r,n))
    Solution1 = (T2 - GG*T1).right_kernel_matrix()[0].list() 
    Solution2 = In_G1*T1 * vector(Solution1)
    V = S(Solution1)
    N = S(Solution2.list())
    ff, re = N.left_quo_rem(-V)
    return ff


def test(totalltests):
    succ = 0
    failure = 0 
    for i in range(totalltests): 
        e = random_small_vector_genenration(m,n,r)    # Generating the error of weight r 
        f = S.random_element(degree=(-1,k-1))  # The message polynomial
        #f = S.random_element(degree = k-1)  # The message polynomial  
        y_list = [f(g[i]) + e[i] for i in range(n)]
        y = vector(y_list)   # 计算被错误干扰的码字
        ff = Decoding_AGCodes(y, r)
        try:
            if (ff == f): 
                succ += 1
            else:
                failure += 1
        except:
            print("Unexpected error", sys.exc_info()[0])
            
    print ("success/totalltests: %d/%d; success rate: %f"  %(succ, totalltests, succ/totalltests))
    print ("failure/totalltests: %d/%d; failure rate: %f" %(failure, totalltests, failure/totalltests))


    
(q,m,n,t,k,r) = (2,45,79,45,5,37)
Fqm.<a> = GF(q**m)
Frob = Fqm.frobenius_endomorphism()
S = OrePolynomialRing(Fqm, Frob, 'x')
# S.<x> = Fqm['x', Frob]

G1 = zero_matrix(Fqm, k+r, k+r)
while G1.is_invertible() == 0:
    g = random_small_vector_genenration(m,n,min(t,m,n)) # Generator of Gabidulin codes
    g_monomials = [g[i]**(q**j) for i in range(n) for j in range(k+r)] 
    G = matrix(Fqm,n,k+r,g_monomials)
    G1 = G.matrix_from_rows(range(k+r))

In_G1 = G1**(-1)
G2 = G.matrix_from_rows(range(k+r,n))
GG = G2 * In_G1


%time test(10)
