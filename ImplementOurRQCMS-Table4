# EG codes + Our RQC-MS (Table 4)

def random_small_vec_gen(n,t):
    B = matrix(Fqm.base_ring(),t,m,0)
    while B.rank() != t:
        B = random_matrix(Fqm.base_ring(),t,m)
    C = matrix(Fqm.base_ring(),n,m,0)
    while C.rank() != t:
        C = random_matrix(Fqm.base_ring(),n,t) * B 
    return vector(Fqm,[C[i] for i in range(n)])

def random_support_gen(t):
    B = matrix(Fqm.base_ring(),t,m,0)
    while B.rank() != t:
        B = matrix(Fqm.base_ring(),[vector(Fqm.random_element()) for i in range(t)])
    return vector(Fqm,[B[i] for i in range(t)])

def random_small_vec_gen_with_support(n,support):
    t = len(list(support))
    C = matrix(Fqm.base_ring(),n,t,0)
    while C.rank() != t:
        C = random_matrix(Fqm.base_ring(),n,t) 
    c = C * support
    return c


def Fold(Codeword, block_length, blocks):
    M = matrix(blocks, block_length, Codeword).transpose()
    return M

def UnFold(FoldCodeword, block_length, blocks):
    return vector(FoldCodeword.transpose().list())

def Vector_Matrix_Dot(Vector, Matrix, block_length, blocks):
    Null_List = []
    for i in range(blocks):
        u = Vector * R(list(Matrix.column(i)))
        Null_List = Null_List + list(u) 
    return Fold(Null_List, block_length, blocks)

def Encoding_EG(Message, SH_Support):
    f = S(Message.list())  # The message polynomial 
    return vector(f.multi_point_evaluation(SH_Support))

def Decoding_EG(Noisy_Word, SH_Support, r): 
    code_length = len(list(SH_Support))
    g_monomials = [SH_Support[i]**(q**j) for i in range(code_length) for j in range(k+r)] 
    SC2 = matrix(Fqm,code_length,k+r,g_monomials) 
    y_monomials = [Noisy_Word[i]**(q**j) for i in range(code_length) for j in range(r+1)] 
    SC1 = matrix(Fqm, code_length, r+1,y_monomials) 
    SC = block_matrix(Fqm, 1, 2, [SC1,SC2])
    #Solution = list(SC.right_kernel().random_element())  
    Solution = SC.right_kernel_matrix()[0].list()  
    V = S(Solution[0:r+1])
    N_vector = vector(Solution[r+1:k+2*r+1])
    N = S(list(N_vector))
    ff,re = N.left_quo_rem(-V)
    return vector(ff.list())
    
# Key Generation
def Blockwise_RQC_MS_KGen(q, m, n, k, w_x, w_y, w_r1, w_r2, w_e):
    h = R.random_element()
    x = R(list(random_small_vec_gen(n, w_x)))
    y = R(list(random_small_vec_gen(n, w_y)))
    s = x + h*y
    pk = [h, s]; sk = [x, y]
    return pk, sk

# Encryption
def Blockwise_RQC_MS_Enc(Public_Key, Message,SH_Support):  
    support1 = random_support_gen(w_r1)
    r1 = random_small_vec_gen_with_support(n,support1)
    r2 = random_small_vec_gen_with_support(n,support1)
    R1 = Fold(list(r1)+list(r2),n, N)
    
    support2 = random_support_gen(w_r2)
    rr1 = random_small_vec_gen_with_support(n,support2)
    rr2 = random_small_vec_gen_with_support(n,support2)
    R2 = Fold(list(rr1)+list(rr2),n, N)
    
    U = R1 + Vector_Matrix_Dot(Public_Key[0], R2, n, N)
    
    support3 = random_support_gen(w_e)
    e1 = random_small_vec_gen_with_support(n,support3)
    e2 = random_small_vec_gen_with_support(n,support3)
    E = Fold(list(e1)+list(e2),n, N)
    
    VV = Vector_Matrix_Dot(Public_Key[1], R2, n, N) + E
    
    V = Fold(Encoding_EG(Message,SH_Support), n, N) + VV 
    ct = [U, V]
    return ct

# Decryption
def Blockwise_RQC_MS_Dec(Private_Key, Ciphertext, SH_Support, r):  
    U = Ciphertext[0]
    V = Ciphertext[1]
    Fold_Noisy_Word = V - Vector_Matrix_Dot(Private_Key[1], U, n, N)
    Noisy_Word = UnFold(Fold_Noisy_Word, n, N)
    return Decoding_EG(Noisy_Word, SH_Support,r)


# EG + RQC-MS
#(q,m,n,k,N,w_x,w_y,w_r1,w_r2,w_e) = (2,47,53,6,2,4,4,4,4,4) # EG + Our RQC-MS -128
#(q,m,n,k,N,w_x,w_y,w_r1,w_r2,w_e) = (2,67,71,4,2,5,5,5,5,6) # EG + Our RQC-MS -192
(q,m,n,k,N,w_x,w_y,w_r1,w_r2,w_e) = (2,73,89,5,2,5,5,6,5,8) # EG + Our RQC-MS -256


Fqm = GF(q**m)
P1  = GF(q)['XX'].irreducible_element(n, algorithm = "minimal_weight")
P.<XX> = Fqm[]
R.<X> = P.quotient(P1)
Frob = Fqm.frobenius_endomorphism()
S = OrePolynomialRing(Fqm, Frob, 'x')

%time Public_Key, Private_Key = Blockwise_RQC_MS_KGen(q, m, n, k, w_x, w_y, w_r1, w_r2, w_e)

Message = random_vector(Fqm, k);  g = random_small_vec_gen(N*n,min(n, m))
%time Ciphertext = Blockwise_RQC_MS_Enc(Public_Key, Message, g)

r = w_x*w_r2 + w_y*w_r1 + w_e
%time Message_test = Blockwise_RQC_MS_Dec(Private_Key, Ciphertext, g, r)

# check correctness
Message_test == Message
