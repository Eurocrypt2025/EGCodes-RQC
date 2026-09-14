############################################################
# Welch-Berlekamp without Loidreau's improved technique proposed in DCC 2018 (section 4.4.2) 
############################################################

def random_small_vec_gen(n,t): 
    B = matrix(Fqm.base_ring(),t,m,0)
    while B.rank() != t:
        B = matrix(Fqm.base_ring(), [vector(Fqm.random_element()) for i in range(t)])
    C = matrix(Fqm.base_ring(),n,m,0)
    while C.rank() != t: 
        C = matrix(Fqm.base_ring(),n,t,[Fqm.base_ring().random_element() for _ in range(n*t)]) * B  
    return vector(Fqm, [C[i] for i in range(n)])

# interpolate:  A(g_i)=0;  I(g_i)=y_i
def interpolate_vect_and_zero(g, y, k):
    A = S(1)
    I = S(0)
    for i in range(k):
        val_A = A(g[i])
        val_I = I(g[i])
        tmp = (y[i] - val_I) / val_A       
        I = I + tmp * A
        A = x * A - (val_A ** (q - 1)) * A   
    return A, I


def reconstruction(g, y, k, n):
    """
    Algorithm 5 in DCC 2018, No using Loidreau's improved technique proposed  (section 4.4.2) 
    Returns (N1, V1, A, I) where (N1, V1) is a solution to LR
    """
    t = (n - k) // 2

    # Build A, I on first k points
    A, I = interpolate_vect_and_zero(g, y, k)

    N0 = A; V0 = S(0)
    N1 = I; V1 = S(1)

    # Discrepancy vectors
    u0 = [A(gi) - V0(yi) for gi, yi in zip(g, y)]
    u1 = [I(gi) - V1(yi) for gi, yi in zip(g, y)]

    # Make mutable copies for swapping
    g_list = list(g)
    y_list = list(y)

    i = k
    while i < n:
        # ----- Secondary loop -----
         # ----- Type 3 handling (u0[i] != 0 and u1[i] == 0) -----
        if  u0[i] != 0 and u1[i] == 0:
            # Find next j > i with u1[j] != 0 or (u0[j]==0 and u1[j]==0)
            j = i + 1
            while j < n and u0[j] != 0 and u1[j] == 0:
                j += 1
            if j == n:
                # No suitable j -> return current (N1, V1) as solution
                return N1, V1
            # Swap positions i and j in all vectors
            g_list[i], g_list[j] = g_list[j], g_list[i]
            y_list[i], y_list[j] = y_list[j], y_list[i]
            u0[i], u0[j] = u0[j], u0[i]
            u1[i], u1[j] = u1[j], u1[i]
            # After swapping, continue with same i (the new u1[i] is non‑zero)
            continue

        # ----- Type 2 handling (u0[i]==0 and u1[i]==0) -----
        # nothing to do for this point (Gabidulin codes do not occur，and AG codes occurs if gi in <g0,g1,...,gk> )
        if u1[i] == 0:
            # newN1 = x * N1
            # newV1 = x * V1
            # newN0 = N0
            # newV0 = V0

            # N0 = newN1
            # V0 = newV1
            # N1 = newN0
            # V1 = newV0

            i += 1
            continue

        # ----- Type 1 handling (u1[i] != 0) -----
        e1 = u1[i] ** (q - 1)      # u1^q / u1
        e2 = u0[i] / u1[i]

        # oldN0, oldV0 = N0, V0

        # # N0' = X*N1 - e1*N1; V0' = X*V1 - e1*V1
        # N0 = x * N1 - e1 * N1
        # V0 = x * V1 - e1 * V1

        # # N1' = N0 - e2*N1; V1' = V0 - e2*V1
        # N1 = oldN0 - e2 * N1
        # V1 = oldV0 - e2 * V1

        # N1' = X*N1 - e1*N1; V1' = X*V1 - e1*V1
        newN1 = x * N1 - e1 * N1
        newV1 = x * V1 - e1 * V1

        # N0' = N0 - e2*N1; V0' = V0 - e2*V1
        newN0 = N0 - e2 * N1
        newV0 = V0 - e2 * V1

        N0 = newN1
        V0 = newV1
        N1 = newN0
        V1 = newV0

        # Update discrepancy vectors for indices > i
        for j in range(i + 1, n):
            new_u0 = (u1[j] ** q) - e1 * u1[j]
            new_u1 = u0[j] - e2 * u1[j]
            u0[j] = new_u0
            u1[j] = new_u1

        i += 1

    return N1, V1

def Welch_Berlekamp_decode(g, y, k, n):
    """
    The decoding of (Augmented) Gabidulin using (Welch-Berlekamp like algorithm)。
    Input:
      g : generator
      y : recieved vector
      k : message length
    Output:
      f : message polynomial (degree <= k-1)
    """
    # reconstruction
    N, W = reconstruction(g, y, k, n)
    # recover f
    Qu, Re = N.left_quo_rem(W)
    return Qu


def test(totalltests):
    succ = 0
    failure = 0
    for npair in range(totalltests):
        e = random_small_vec_gen(n, r)
        y = c + e
        f_rec = Welch_Berlekamp_decode(g, y, k, n)
        try:
            if (f == f_rec): 
                succ += 1
            else:
                failure += 1
        except:
            print("Unexpected error", sys.exc_info()[0])
            
    print ("success/totalltests: %d/%d; success rate: %f" % (succ,totalltests,succ/totalltests))
    print ("failure/totalltests: %d/%d; failure rate: %f" % (failure,totalltests,failure/totalltests))



# Compute Theoretical and Simulated DFR by Welch_Berlekamp like algorithm
# increase m
#(q,m,n,t,k,r) = (2,31,41,31,9,16)   # DFR: 2**(-5)
#(q,m,n,t,k,r) = (2,32,41,32,9,16)   # DFR: 2**(-6)
#(q,m,n,t,k,r) = (2,33,41,33,9,16)   # DFR: 2**(-7)
#(q,m,n,t,k,r) = (2,34,41,34,9,16)   # DFR: 2**(-8)
#(q,m,n,t,k,r) = (2,35,41,35,9,16)   # DFR: 2**(-9)
# TheoreticalDFR = [0.0313, 0.0156, 0.0078, 0.0039, 0.0020]
# SimulatedDFR = [0.0150, 0.0076, 0.0038, 0.0017, 0.0008]

# increase t
#(q,m,n,t,k,r) = (2,35,41,30,9,16)   # DFR: 2**(-4)
#(q,m,n,t,k,r) = (2,35,41,31,9,16)   # DFR: 2**(-5)
#(q,m,n,t,k,r) = (2,35,41,32,9,16)   # DFR: 2**(-6)
#(q,m,n,t,k,r) = (2,35,41,33,9,16)   # DFR: 2**(-7)
#(q,m,n,t,k,r) = (2,35,41,34,9,16)   # DFR: 2**(-8)
# TheoreticalDFR = [0.0625, 0.0313, 0.0156, 0.0078, 0.0039]
# SimulatedDFR = [0.0310, 0.0150, 0.0084, 0.0036, 0.0018] 

# increase t < n < m; increase t
#(q,m,n,t,k,r) = (2,29,26,16,5,10)   # DFR: 2**(-2)
#(q,m,n,t,k,r) = (2,29,26,17,5,10)   # DFR: 2**(-4)
#(q,m,n,t,k,r) = (2,29,26,18,5,10)   # DFR: 2**(-6)
#(q,m,n,t,k,r) = (2,29,26,19,5,10)   # DFR: 2**(-8)
#(q,m,n,t,k,r) = (2,29,26,20,5,10)   # DFR: 2**(-10)
# TheoreticalDFR = [0.2500, 0.0625, 0.0156, 0.0039, 0.00098]
# SimulatedDFR = [0.1320, 0.0372, 0.0098, 0.0025, 0.00058] 

# Decoding up to the RGV bound (k > r)
#(q,m,n,t,k,r) = (2,30,37,30,23,7)   # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,30,38,30,22,8)   # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,30,39,30,21,9)   # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,30,40,30,20,10)  # DFR: 2**(1) 
#(q,m,n,t,k,r) = (2,30,41,30,19,11)  # DFR: 2**(1)
# TheoreticalDFR = [1, 1, 1, 1, 1]
# SimulatedDFR = [0.7120, 0.7087, 0.7100, 0.7105, 0.7141] 

# Decoding up to the RGV bound (k < r)
#(q,m,n,t,k,r) = (2,21,34,21,8,13)  # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,22,36,22,8,14)  # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,23,38,23,8,15)  # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,24,40,24,8,16)  # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,25,42,25,8,17)  # DFR: 2**(1)
# TheoreticalDFR = [1, 1, 1, 1, 1]
# SimulatedDFR = [0.7124, 0.7120, 0.7118, 0.7120, 0.7125]


# increase n
#(q, m, n, t, k, r) = (2, 27, 41, 27, 9, 16)  
#(q, m, n, t, k, r) = (2, 27, 42, 27, 9, 16) 
(q, m, n, t, k, r) = (2, 27, 43, 27, 9, 16) 
#(q, m, n, t, k, r) = (2, 27, 44, 27, 9, 16)
#(q, m, n, t, k, r) = (2, 27, 45, 27, 9, 16)
# TheoreticalDFR = [0.5000, 0.0625, 0.0078, 0.00098, 0.00012]
# SimulatedDFR = [0.2290, 0.0384, 0.0057, 0.00074, 0.00003] 



Fqm.<a> = GF(q**m)
Frob = Fqm.frobenius_endomorphism()
# S = OrePolynomialRing(Fqm, Frob, 'x')
S.<x> = Fqm['x', Frob]

g = random_small_vec_gen(n, min(m, n, t))

#f = (a*x^0 + a^3*x + a^7*x^2)
Message = random_vector(Fqm, k) 
f = S(list(Message))
c = vector(f.multi_point_evaluation(g))

%time test(1000)
