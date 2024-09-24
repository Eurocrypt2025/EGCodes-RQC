# Decoding any syndrome by EG codes to errors of weight the RGV bound 

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
        s = random_vector(Fqm,n-k)  # random syndrome
        y = H.solve_right(s)       
        y_monomials = [y[i]**(q**j) for i in range(n) for j in range(r+1)] 
        A1 = matrix(Fqm,n,r+1,y_monomials) 
        A = block_matrix(Fqm,1,2,[A1,A2])
        Solution = A.right_kernel_matrix()[0].list()  # Decoding step
        V = S(Solution[0: r+1])
        N_vector = vector(Solution[r+1: k+2*r+1])
        N = S(list(N_vector))
        ff, re = N.left_quo_rem(-V)
        e = y - vector(ff.multi_point_evaluation(g))
        try:
            if (rank_weight(e, Fqm.base_ring()) <=  r) and (H*y == H*e == s): 
                succ += 1
            else:
                failure += 1
        except:
            print("Unexpected error", sys.exc_info()[0])
            
    print ("success/totalltests: %d/%d; success rate: %f" % (succ,totalltests,succ/totalltests))
    print ("failure/totalltests: %d/%d; failure rate: %f" % (failure,totalltests,failure/totalltests))
    
    
# Decoding up to the RGV bound (k > r)
(q,m,n,t,k,r) = (2,30,37,30,23,7)   # Can decode weight r = 8
#(q,m,n,t,k,r) = (2,30,38,30,22,8)   # Can decode weight r = 9
#(q,m,n,t,k,r) = (2,30,39,30,21,9)   # Can decode weight r = 10
#(q,m,n,t,k,r) = (2,30,40,30,20,10)  # Can decode weight r = 11
#(q,m,n,t,k,r) = (2,30,41,30,19,11)  # Can decode weight r = 12
# RGVBound = [7, 8, 9, 10, 11]
# RSBound = [12, 13, 14, 16, 17]


# Decoding up to the RGV bound (k < r)
#(q,m,n,t,k,r) = (2,21,34,21,8,13)  
#(q,m,n,t,k,r) = (2,22,36,22,8,14)  
#(q,m,n,t,k,r) = (2,23,38,23,8,15)  
#(q,m,n,t,k,r) = (2,24,40,24,8,16)  
#(q,m,n,t,k,r) = (2,25,42,25,8,17)  
# RGVBound = [13, 14, 15, 16, 17]
# RSBound = [17, 18, 19, 20, 21]


Fqm.<a> = GF(q**m)
Frob = Fqm.frobenius_endomorphism()
S = OrePolynomialRing(Fqm, Frob, 'x')
# S.<x> = Fqm['x', Frob]
#print(S)

g = random_small_vector_genenration(m, n,min(m,n,t))
g_monomials = [g[i]**(q**j) for i in range(n) for j in range(k+r)] 
A2 = matrix(Fqm, n, k+r, g_monomials) 

G = matrix(Fqm,k,n,[g[i]**(q**j) for j in range(k) for i in range(n)]) # Generator matrix
H = G.right_kernel_matrix()  # Parity-check matrix


%time test(100000)
