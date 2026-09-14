################### The MM Modeling for the Blockwise RD problem ##########################
################################### Two Blocks #########################################
###########################################################################################

from itertools import combinations, combinations_with_replacement

def support_basis_generation(q,m,n,k,r,n1,n2,r1,r2):
    S = matrix(Fqm.base_ring(), r, m, 0)
    while S.rank() != r:
        S = random_matrix(Fqm.base_ring(), r, m)
    return vector(Fqm, [S[i] for i in range(r)])  # S[i]: i-th row of S

def coefficient_matrix_generation(q,m,n,k,r,n1,n2,r1,r2):
    C1 = block_matrix(1, 2, [identity_matrix(r1), random_matrix(Fqm.base_ring(), r1, n1-r1)])
    C2 = block_matrix(1, 2, [identity_matrix(r2), random_matrix(Fqm.base_ring(), r2, n2-r2)])
    C = block_matrix(2, 2, [C1,0,0,C2])
    return C

#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,1,4,6,6,2,2) 
# SME & SNE & SLIE:  210 & 185 & 160; 0 < k+1 <= r1
#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,2,4,6,6,2,2) 
# SME & SNE & SLIE:  126 & 120 & 114; r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,3,4,6,6,2,2) 
# SME & SNE & SLIE:  70 & 70 & 70; r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,4,4,6,6,2,2) 
# SME & SNE & SLIE:  35 & 35 & 35; r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,5,4,6,6,2,2) 
# SME & SNE & SLIE:  15 & 15 & 15; r1 < k+1 <= n1

#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,3,4,8,4,2,2) 
# SME & SNE & SLIE:  70 & 70 & 53; 0 < k+1 <= r1
#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,4,4,8,4,2,2) 
# SME & SNE & SLIE:  35 & 35 & 31; r1 < k+1 <= n1
(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,5,4,8,4,2,2) 
# SME & SNE & SLIE:  15 & 15 & 15; r1 < k+1 <= n1


TME = binomial(n-k-1, r); Unknowns_number = binomial(n1, r1)*binomial(n2, r2)
print("Theoretical Maximal Equations (Fqm, TME) : ", TME)
print("The number of Unknowns: ", Unknowns_number)


Fqm = GF(q**m)

RD_generator_matrix = random_matrix(Fqm, k, n)
while RD_generator_matrix.rank() != k:
    RD_generator_matrix = random_matrix(Fqm, k, n)
RD_parity_check_matrix = RD_generator_matrix.right_kernel_matrix()
message = random_vector(Fqm, k)

# Generate Blockwise Error with Three Blocks
support_basis = support_basis_generation(q,m,n,k,r,n1,n2,r1,r2)  # r  
coefficient_matrix = coefficient_matrix_generation(q,m,n,k,r,n1,n2,r1,r2)  # r * n  matrix
error = support_basis * coefficient_matrix

# RD instance
y = message * RD_generator_matrix + error

# Reduce the RD instance to finding small codeword in the extended code 
generator_matrix = block_matrix(2, 1, [matrix(y), RD_generator_matrix])

parity_check_matrix = generator_matrix.right_kernel_matrix()  # systematic parity check matrix

# U = random_matrix(Fqm, n-k-1, n-k-1)  
# while U.rank() != n-k-1:
#       U = random_matrix(Fqm, n-k-1, n-k-1)
# parity_check_matrix = U*parity_check_matrix   # general parity check matrix

CH_T = coefficient_matrix * parity_check_matrix.transpose()
CH_T_left_kernel = CH_T.left_kernel()
print("Dimension of Left Kernel of CH^T:", CH_T_left_kernel.dimension())


var_names = ["m"+str(i) for i in range(Unknowns_number)]
PR = PolynomialRing(Fqm.base_ring(), var_names)

L1 = list(combinations(list(range(n1)), r1))
L2 = list(combinations(list(range(n1, n)), r2))
L = list(itertools.product(L1, L2))

V = []
for i in range(Unknowns_number):
    LL = L[i]; v = LL[0] + LL[1]
    V.append(v)


# Dictionary for m_variables
m_vars = {}; ctr = 0
for comb in V:
    m_vars[comb] = PR.gens()[ctr]  
    ctr += 1

linear_eqns = []
for J in combinations(list(range(n-k-1)), r):
    sub_parity_check_matrix_T = parity_check_matrix.transpose().matrix_from_columns(J)
    eqn = 0
    for T in V:
        H_T_T = sub_parity_check_matrix_T.matrix_from_rows(T).det()
        eqn += m_vars[T] * H_T_T
    linear_eqns.append(eqn)

SME = len(linear_eqns); SNE = len(linear_eqns) - linear_eqns.count(0)
print("Simulated Maximal Equations (Fqm, SME): ", SME)
print("Simulated Non-zero Equations (Fqm, SNE): ", SNE)

# Compute the ranks of the linear equations system
def computeRanks(linear_eqns):
    s = Sequence(linear_eqns)
    if len(linear_eqns) > 0:
        M, Mon = s.coefficient_matrix()  # Coefficients matrix and unknowns for system
        rows = M.nrows(); cols = M.ncols(); Rank = M.rank()
        return rows, cols, Rank
    else:
        print("no equations")
    return 0, 0, 0

print()
print("--- RANK EXPERIMENTS ---")
rows, cols, Rank = computeRanks(linear_eqns) 
print("rows : %d, cols : %d, rank : %d" %(rows, cols, Rank))
print("SME & SNE & SLIE: ", rows, "&", SNE, "&", Rank)
