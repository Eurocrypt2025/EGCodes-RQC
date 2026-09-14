# Test decoding failure rate of AG codes
# Decoding errors of weight exactly error_weight w < r
# Standard Gaussian elimination 

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
        e = random_small_vector_genenration(m, n, w)    # Generating the error e of weight w
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
# r <= min(n-k/2,m-k)


# Compute Theoretical DFR by Theorem 3 for code parameters in Table 5
# Decrease w (t = m < n)
#(q,m,n,t,k,r,w) = (2, 27, 41, 27, 9, 16, 16)   # DFR: 2**(-1)
#(q,m,n,t,k,r,w) = (2, 27, 41, 27, 9, 16, 15)  # DFR: 2**(-4)
#(q,m,n,t,k,r,w) = (2, 27, 41, 27, 9, 16, 14)  # DFR: 2**(-7)
#(q,m,n,t,k,r,w) = (2, 27, 41, 27, 9, 16, 13)  # DFR: 2**(-10)
#(q,m,n,t,k,r,w) = (2, 27, 41, 27, 9, 16, 12)  # DFR: 2**(-13)
# TheoreticalDFR = [0.5000, 0.0625, 0.0078, 0.00098, 0.00012]
# SimulatedDFR = [0.2312, 0.0380, 0.0050,  0.00062,  0.00007] 

# Decrease w, (t < n < m)
#(q,m,n,t,k,r,w) = (2, 29, 26, 16, 5, 10, 10)   # DFR: 2**(-2)
#(q,m,n,t,k,r,w) = (2, 29, 26, 16, 5, 10, 9)  # DFR: 2**(-4)
#(q,m,n,t,k,r,w) = (2, 29, 26, 16, 5, 10, 8)  # DFR: 2**(-6)
#(q,m,n,t,k,r,w) = (2, 29, 26, 16, 5, 10, 7)   # DFR: 2**(-8)
#(q,m,n,t,k,r,w) = (2, 29, 26, 16, 5, 10, 6)  # DFR: 2**(-10)
# TheoreticalDFR = [0.2500, 0.0625, 0.0156, 0.0039, 0.00098]
# SimulatedDFR = [0.1326, 0.0362, 0.0097, 0.0026, 0.00069] 

# Decrease w, Decoding up to the RGV bound (k > r)
# Hash-Sign 
#(q,m,n,t,k,r,w) = (2, 30, 37, 30, 23, 7, 7)   # DFR: 2**(1)
#(q,m,n,t,k,r,w) = (2, 30, 37, 30, 23, 7, 6)   # DFR: 2**(0)
#(q,m,n,t,k,r,w) = (2, 30, 37, 30, 23, 7, 5)   # DFR: 2**(-1)
#(q,m,n,t,k,r,w) = (2, 30, 37, 30, 23, 7, 4)   # DFR: 2**(-2)
#(q,m,n,t,k,r,w) = (2, 30, 37, 30, 23, 7, 3)  # DFR: 2**(-3)
# TheoreticalDFR = [1,  1,  0.5,  0.2500, 0.125]
# SimulatedDFR = [0.7122, 0.4197, 0.2217, 0.1145, 0.0549] 

# Decrease w,  Decoding up to the RGV bound (k < r)
#(q,m,n,t,k,r,w) = (2, 21, 34, 21, 8, 13, 13)  # DFR: 2**(1)
#(q,m,n,t,k,r,w) = (2, 21, 34, 21, 8, 13, 12)  # DFR: 2**(0)
#(q,m,n,t,k,r,w) = (2, 21, 34, 21, 8, 13, 11) # DFR: 2**(-1)
#(q,m,n,t,k,r,w) = (2, 21, 34, 21, 8, 13, 10)  # DFR: 2**(-2)
#(q,m,n,t,k,r,w) = (2, 21, 34, 21, 8, 13, 9)  # DFR: 2**(-3)
# TheoreticalDFR = [1, 1, 0.5, 0.2500, 0.125]
# SimulatedDFR = [0.7125, 0.4234, 0.2312, 0.1203, 0.0608]


Fqm.<a> = GF(q**m)
Frob = Fqm.frobenius_endomorphism()
S = OrePolynomialRing(Fqm, Frob, 'x')
# S.<x> = Fqm['x', Frob]
#print(S)
     
g = random_small_vector_genenration(m,n,min(t,m,n))  # Generator of SH codes
g_monomials = [g[i]**(q**j) for i in range(n) for j in range(k+r)] 
A2 = matrix(Fqm, n, k+r, g_monomials) 

%time test(100000)
