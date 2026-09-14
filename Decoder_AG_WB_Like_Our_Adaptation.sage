# ============================================================
# Welch--Berlekamp-like reconstruction for AG codes
# Pair (N,W) represents N(g_i) = W(y_i)
# In the paper notation: N = u,  W = v.
# The shifted degree is max(deg_q N, deg_q W + k - 1).
# ============================================================


def _ore_degree(P):
    """
    q-degree of an Ore polynomial. We use a sufficiently small value for the zero polynomial.
    """
    if P == 0:
        return -10**9
    return ZZ(P.degree())

def WB_shifted_degree(pair, k):
    """
    (0,k-1)-shifted degree: sdeg(N,W) = max(deg N, deg W + k - 1).
    """
    N, W = pair
    return max(_ore_degree(N), _ore_degree(W) + k - 1)

def WB_leading_position(pair, k):
    """
    Shifted leading position.
    0 : N-component
    1 : W-component
    Ties are resolved in favor of the W-component.
    """
    N, W = pair
    dN = _ore_degree(N)
    dW = _ore_degree(W) + k - 1
    if dW >= dN:
        return 1
    else:
        return 0

def WB_leading_key(pair, k):
    """
    Key used to compare shifted leading monomials.
    Smaller key = smaller shifted leading monomial.
    """
    return (WB_shifted_degree(pair, k), WB_leading_position(pair, k))

def WB_discrepancy(pair, gi, yi):
    """
    Discrepancy: delta = N(g_i) - W(y_i).
    """
    N, W = pair
    return N(gi) - W(yi)

def WB_annihilator_one_element(delta):
    """
    Degree-one q-annihilator of <delta>_Fq.
    In the Ore representation, A_delta = x - theta(delta)/delta,
    since A_delta(delta) = theta(delta) - theta(delta)/delta * delta = 0.
    """
    if delta == 0:
        raise ValueError("delta must be nonzero")
    gamma = Frob(delta) / delta
    return S.gen() - gamma

def WB_AG_LRReconstruction(g, y, k, debug=False):
    """
    Welch--Berlekamp-like reconstruction for AG codes.
    INPUT:
        g -- evaluation vector
        y -- received word
        k -- code dimension

    OUTPUT:
        (N,W) such that N(g_i) = W(y_i) for every i, and (N,W) has minimum shifted degree among the interpolation-module elements.

    No F_q-linear independence of all coordinates g_i is required.
    """

    n = len(g)

    if len(y) != n:
        raise ValueError("g and y must have the same length")

    # --------------------------------------------------------
    # Initial basis of M_0 = L[x]^2.
    # S.one() is the identity q-polynomial in the Ore
    # representation.
    # --------------------------------------------------------

    b = [(S.one(), S.zero()),   # b_0 = (1,0)
         (S.zero(), S.one())    # b_1 = (0,1)
    ]

    # --------------------------------------------------------
    # Process interpolation constraints one by one: N(g_i) = W(y_i).
    # --------------------------------------------------------
    for i in range(n):
        delta = [WB_discrepancy(b[0], g[i], y[i]), WB_discrepancy(b[1], g[i], y[i])]

        # Both basis elements already satisfy the new constraint.
        if delta[0] == 0 and delta[1] == 0:
            if debug:
                print("i =", i,
                      "type = redundant",
                      "sdeg =", [WB_shifted_degree(B, k) for B in b])
            continue

        # ----------------------------------------------------
        # Pivot selection:
        # among nonzero-discrepancy rows, choose the one with smaller shifted leading monomial.
        # ----------------------------------------------------
        candidates = [j for j in range(2) if delta[j] != 0]
        p = min(candidates, key=lambda j: WB_leading_key(b[j], k))
        s = 1 - p

        bp_old = b[p]
        bs_old = b[s]

        dp = delta[p]
        ds = delta[s]

        # ----------------------------------------------------
        # Non-pivot update: b_s <- b_s - (delta_s/delta_p) b_p.
        # This kills the new discrepancy.
        # ----------------------------------------------------
        beta = ds / dp
        b_s_new = (bs_old[0] - beta * bp_old[0], bs_old[1] - beta * bp_old[1] )

        # ----------------------------------------------------
        # Pivot update:   b_p <- A_delta_p o b_p
        # where A_delta_p(delta_p)=0.
        # Ore multiplication on the left corresponds
        # to symbolic composition on the left.
        # ----------------------------------------------------
        A_delta = WB_annihilator_one_element(dp)
        b_p_new = (A_delta * bp_old[0], A_delta * bp_old[1])
        b[p] = b_p_new
        b[s] = b_s_new

        if debug:
            # Check that both basis rows satisfy all processed
            # interpolation constraints.
            for row in range(2):
                Nrow, Wrow = b[row]

                assert all(Nrow(g[j]) == Wrow(y[j]) for j in range(i + 1))

            print("i =", i,
                  "pivot =", p,
                  "sdeg =", [WB_shifted_degree(B, k) for B in b],
                  "lp =", [WB_leading_position(B, k) for B in b])

    # --------------------------------------------------------
    # Return the minimum shifted-degree basis element.
    # --------------------------------------------------------
    j_min = min(range(2), key=lambda j: WB_leading_key(b[j], k))

    return b[j_min], b

def RankWeight_qPolynomial(v):
    """
    Rank weight over the base field associated with the Frobenius automorphism used by S.
    """
    vv = list(v)
    if all(a == 0 for a in vv):
        return 0
    A = S.minimal_vanishing_polynomial(vv)
    return ZZ(A.degree())

def WelchBerlekamp_AG_Decoding(n, k, r, g, y, debug=False, return_details=False):
    """
    Welch--Berlekamp-like decoder for AG codes.
    Returns:   f  if decoding succeeds;
               None  if decoding fails.
    If return_details=True, returns (f, e, N, W, basis).
    """

    # --------------------------------------------------------
    # Step 1: WB-like LR reconstruction
    # --------------------------------------------------------
    (N, W), basis = WB_AG_LRReconstruction(g, y, k, debug=debug)

    # --------------------------------------------------------
    # Step 2: W must be nonzero.
    # Under k+r <= t and the good-kernel condition, this cannot happen, but we explicitly test it.
    # --------------------------------------------------------
    if W == 0:
        if debug:
            print("Decoding failure: W = 0")
        return None

    # --------------------------------------------------------
    # Step 3: symbolic left division   N = W o f + R.
    # --------------------------------------------------------
    f, rem = N.left_quo_rem(W)
    if rem != 0:
        if debug:
            print("Decoding failure: nonzero division remainder")
        return None

    # --------------------------------------------------------
    # Step 4: message-degree check
    # --------------------------------------------------------
    if f != 0 and f.degree() > k - 1:
        if debug:
            print("Decoding failure: deg(f) > k-1")
        return None

    # --------------------------------------------------------
    # Step 5: reconstruct the codeword and error
    # --------------------------------------------------------
    codeword = vector(Fqm, f.multi_point_evaluation(list(g)))
    e = vector(Fqm, y) - codeword; wt = RankWeight_qPolynomial(e)

    # --------------------------------------------------------
    # Step 6: validation
    # --------------------------------------------------------
    if wt > r:
        if debug:
            print("Decoding failure: rank weight =", wt, "> r =", r)
        return None

    if return_details:
        return f, e, N, W, basis

    return f

def Encoding_AGabidulin(Message, AG_Support):
    f = S(Message.list())             #  The message polynomial 
    return vector(f.multi_point_evaluation(AG_Support))


def random_small_vector_genenration(Extension, Length, Weight):
    B = matrix(Fqm.base_ring(), Weight, Extension, 0)
    while B.rank() != Weight:
        B = random_matrix(Fqm.base_ring(),Weight, Extension)
    C = matrix(Fqm.base_ring(), Length, Extension,0)
    while C.rank() != Weight:
        C = random_matrix(Fqm.base_ring(), Length, Weight) * B
    return vector(Fqm,[C[i] for i in range(Length)])



def test(total_tests):
    succ = 0
    failure = 0
    errors = 0
    for _ in range(total_tests):
        e = random_small_vector_genenration(m, n, r)
        y = Codeword + e
        try:
            result = WelchBerlekamp_AG_Decoding(n, k, r, g, y, debug=False,return_details=True)

            if result is None:
                failure += 1
                continue

            f_decoded, e_decoded, _, _, _ = result
            if (vector(f_decoded.padded_list(k)) == Message and e_decoded == e):
                succ += 1
            else:
                failure += 1

        except Exception as err:
            errors += 1
            print(f"Unexpected error: {type(err).__name__}: {err}")

    tested = succ + failure

    print(f"successful decodings : {succ}")
    print(f"decoding failures    : {failure}")
    print(f"program errors        : {errors}")

    if tested > 0:
        success_rate = float(succ) / float(tested)
        failure_rate = float(failure) / float(tested)
        print(f"success rate          : {success_rate:.6f}")
        print(f"failure rate (DFR)    : {failure_rate:.6f}")


# Compute Theoretical DFR (Theorem 3) and Simulated DFR for code parameters in Table 5
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
(q,m,n,t,k,r) = (2,29,26,20,5,10)   # DFR: 2**(-10)
# TheoreticalDFR = [0.2500, 0.0625, 0.0156, 0.0039, 0.00098]
# SimulatedDFR = [0.1320, 0.0372, 0.0098, 0.0025, 0.00058] 

# Decoding up to the RGV bound (k > r)
# Hash-Sign 
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


Fqm.<a> = GF(q**m)
Frob = Fqm.frobenius_endomorphism()
# S = OrePolynomialRing(Fqm, Frob, 'x')
S.<x> = Fqm['x', Frob]

Message = random_vector(Fqm, k) 
g = random_small_vector_genenration(m, n, min(m, n, t))
# g1 = random_small_vec_gen(t, t);  g2 = zero_vector(Fqm, n-t);  g = vector(g1.list() + g2.list())
Codeword = Encoding_Gabidulin(Message, g)

%time test(100000)
