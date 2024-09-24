# Test decoding failure probability of EG codes
# Decoding errors of weight exactly r


from sage.coding.linear_rank_metric import rank_weight

def random_small_vector_genenration(Extension, Length, Weight):
    B = matrix(Fqm.base_ring(), Weight, Extension, 0)
    while B.rank() != Weight:
        B = random_matrix(Fqm.base_ring(),Weight, Extension)
    C = matrix(Fqm.base_ring(), Length, Extension,0)
    while C.rank() != Weight:
        C = random_matrix(Fqm.base_ring(), Length, Weight) * B
    return vector(Fqm,[C[i] for i in range(Length)])


def test(totalltests):
    succ = 0
    failure = 0
    for npair in range(totalltests):
        f = S.random_element(degree=(-1,k-1))  # The message polynomial
        e = random_small_vector_genenration(m, n, r)    # Generating the error e of weight r 
        y = vector([f(g[i]) + e[i] for i in range(n)])    
        y_monomials = [y[i]**(q**j) for i in range(n) for j in range(r+1)] 
        A1 = matrix(Fqm, n, r+1, y_monomials) 
        A = block_matrix(Fqm, 1, 2, [A1,A2])
        Solution = A.right_kernel_matrix()[0].list()  
        #Solution = list(A.right_kernel().random_element()) 
        V = S(Solution[0: r+1])
        N_vector = vector(Solution[r+1: k+2*r+1])
        N = S(list(N_vector))
        ff, re = N.left_quo_rem(-V)
        try:
            if (ff ==  f) and (y - vector(ff.multi_point_evaluation(g)) == e): 
                succ += 1
            else:
                failure += 1
        except:
            print("Unexpected error", sys.exc_info()[0])
            
    print ("success/totalltests: %d/%d; success rate: %f" % (succ,totalltests,succ/totalltests))
    print ("failure/totalltests: %d/%d; failure rate: %f" % (failure,totalltests,failure/totalltests))


# t = |g|
# t >= k + r  
# n >= k + 2r 
# Ensure the Rank of n * (k + 2r + 1) to be  k+2r 
#  r <= min(n-k/2,m-k)

# Test DFR for code parameters in Table 2
# increase q
#(q,m,n,t,k,r) = (2,5,7,5,2,2)    # DFR: 2**(-2.0)
#(q,m,n,t,k,r) = (3,5,7,5,2,2)    # DFR: 2**(-5.3)
#(q,m,n,t,k,r) = (5,5,7,5,2,2)    # DFR: 2**(-8.3)
#(q,m,n,t,k,r) = (7,5,7,5,2,2)    # DFR: 2**(-10.2)
#(q,m,n,t,k,r) = (11,5,7,5,2,2)   # DFR: 2**(-12.8)
# TheoreticalDFR = [0.2500, 0.0247, 0.0032, 0.00085, 0.00014]
# SimulatedDFR = [0.0579, 0.0123, 0.0016, 0.00034, 0.00004]


# increase m
#(q,m,n,t,k,r) = (2,31,41,31,9,16)   # DFR: 2**(-5)
#(q,m,n,t,k,r) = (2,32,41,32,9,16)   # DFR: 2**(-6)
#(q,m,n,t,k,r) = (2,33,41,33,9,16)   # DFR: 2**(-7)
#(q,m,n,t,k,r) = (2,34,41,34,9,16)   # DFR: 2**(-8)
#(q,m,n,t,k,r) = (2,35,41,35,9,16)   # DFR: 2**(-9)
# TheoreticalDFR = [0.0313, 0.0156, 0.0078, 0.0039, 0.0020]
# SimulatedDFR = [0.0154, 0.0078, 0.0037, 0.0018, 0.0009]

# increase n
#(q,m,n,t,k,r) = (2,27,41,27,9,16)   # DFR: 2**(-1)
#(q,m,n,t,k,r) = (2,27,42,27,9,16)   # DFR: 2**(-4)
#(q,m,n,t,k,r) = (2,27,43,27,9,16)   # DFR: 2**(-7)
#(q,m,n,t,k,r) = (2,27,44,27,9,16)   # DFR: 2**(-10)
#(q,m,n,t,k,r) = (2,27,45,27,9,16)   # DFR: 2**(-13)
# TheoreticalDFR = [0.5000, 0.0625, 0.0078, 0.00098, 0.00012]
# SimulatedDFR = [0.2312, 0.0388, 0.0052, 0.00085, 0.00008] 

# increase t
#(q,m,n,t,k,r) = (2,35,41,30,9,16)   # DFR: 2**(-4)
#(q,m,n,t,k,r) = (2,35,41,31,9,16)   # DFR: 2**(-5)
#(q,m,n,t,k,r) = (2,35,41,32,9,16)   # DFR: 2**(-6)
#(q,m,n,t,k,r) = (2,35,41,33,9,16)   # DFR: 2**(-7)
#(q,m,n,t,k,r) = (2,35,41,34,9,16)   # DFR: 2**(-8)
# TheoreticalDFR = [0.0625, 0.0313, 0.0156, 0.0078, 0.0039]
# SimulatedDFR = [0.0300, 0.0151, 0.0083, 0.0036, 0.0018] 

# increase t < n < m; increase t
#(q,m,n,t,k,r) = (2,29,26,16,5,10)   # DFR: 2**(-2)
#(q,m,n,t,k,r) = (2,29,26,17,5,10)   # DFR: 2**(-4)
#(q,m,n,t,k,r) = (2,29,26,18,5,10)   # DFR: 2**(-6)
#(q,m,n,t,k,r) = (2,29,26,19,5,10)   # DFR: 2**(-8)
#(q,m,n,t,k,r) = (2,29,26,20,5,10)   # DFR: 2**(-10)
# TheoreticalDFR = [0.2500, 0.0625, 0.0156, 0.0039, 0.00098]
# SimulatedDFR = [0.1326, 0.0371, 0.0098, 0.0024, 0.00058] 

# Decoding up to the RGV bound (k > r)
# Hash-Sign 
#(q,m,n,t,k,r) = (2,30,37,30,23,7)   # DFR: 2**(0)
#(q,m,n,t,k,r) = (2,30,38,30,22,8)   # DFR: 2**(0)
#(q,m,n,t,k,r) = (2,30,39,30,21,9)   # DFR: 2**(0)
#(q,m,n,t,k,r) = (2,30,40,30,20,10)  # DFR: 2**(0)   
#(q,m,n,t,k,r) = (2,30,41,30,19,11)  # DFR: 2**(0)
# TheoreticalDFR = [1, 1, 1, 1, 1]
# SimulatedDFR = [0.7122, 0.7086, 0.7102, 0.7109, 0.7142] 

# Decoding up to the RGV bound (k < r)
#(q,m,n,t,k,r) = (2,21,34,21,8,13)  # DFR: 2**(0)
#(q,m,n,t,k,r) = (2,22,36,22,8,14)  # DFR: 2**(0)
#(q,m,n,t,k,r) = (2,23,38,23,8,15)  # DFR: 2**(0)
#(q,m,n,t,k,r) = (2,24,40,24,8,16)  # DFR: 2**(0)
#(q,m,n,t,k,r) = (2,25,42,25,8,17)  # DFR: 2**(0)
# TheoreticalDFR = [1, 1, 1, 1, 1]
# SimulatedDFR = [0.7125, 0.7119, 0.7118, 0.7122, 0.7123]

# Gabidulin codes: t = n <= m
#(q,m,n,t,k,r) = (2,27,27,27,7,10)   # DFR: 0
#(q,m,n,t,k,r) = (2,28,27,27,7,10)   # DFR: 0
#(q,m,n,t,k,r) = (2,29,27,27,7,10)   # DFR: 0
#(q,m,n,t,k,r) = (2,30,27,27,7,10)   # DFR: 0
#(q,m,n,t,k,r) = (2,31,27,27,7,10)   # DFR: 0


Fqm.<a> = GF(q**m)
Frob = Fqm.frobenius_endomorphism()
S = OrePolynomialRing(Fqm, Frob, 'x')
# S.<x> = Fqm['x', Frob]
#print(S)
     
g = random_small_vector_genenration(m,n,min(t,m,n))  # Generator of SH codes
g_monomials = [g[i]**(q**j) for i in range(n) for j in range(k+r)] 
A2 = matrix(Fqm, n, k+r, g_monomials) 


%time test(100000)
