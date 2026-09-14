# Implement Muti-UR-AG 

def random_small_space_gen(t):
    B = matrix(Fqm.base_ring(),t,m,0)
    while B.rank() != t:
        B = matrix(Fqm.base_ring(),[vector(Fqm.random_element()) for i in range(t)])
    return B.row_space()

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


def Fold(Codeword, blocks, block_length): # blocks * block_length
    M = matrix(blocks, block_length, Codeword)
    return M

def UnFold(FoldCodeword):
    return vector(FoldCodeword.list())

def homogenous_matrix_gen(k,n,d): # generate k * n homogenous matrix of dimension d
    F = random_small_space_gen(d)
    H_list = [Fqm(list(F.random_element())) for _ in range(n*k)]
    return matrix(Fqm,k,n,H_list)

def homogenous_matrix_gen_with_support(n1,n2,support):
    t = len(list(support))
    B = matrix(Fqm.base_ring(),[vector(support[i]) for i in range(t)])
    F = B.row_space()
    H_list = [Fqm(list(F.random_element())) for _ in range(n1*n2)]
    return matrix(Fqm,n1,n2,H_list)

def homogenous_matrix_gen_with_support_one(n1,n2,support):
    t = len(list(support))
    support = support[0]**(-1) * support
    B = matrix(Fqm.base_ring(),[vector(support[i]) for i in range(t)])
    F = B.row_space()
    H_list = [Fqm(list(F.random_element())) for _ in range(n1*n2)]
    return matrix(Fqm,n1,n2,H_list)

def Encoding_AG(Message, SH_Support):
    f = S(Message.list())  # The message polynomial 
    return vector(f.multi_point_evaluation(SH_Support))

def Decoding_AG(Noisy_Word, SH_Support, r): 
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
def MutiURAG_KGen(q,m,n,N1,N2,k,t,w_1,w_2):
    H = random_matrix(Fqm, n, n)
    Support = random_support_gen(w_1)
    X = homogenous_matrix_gen_with_support_one(n,N1,Support)
    Y = homogenous_matrix_gen_with_support_one(n,N1,Support)
    S = X + H*Y
    pk = [H, S]; sk = [X, Y]
    return pk, sk

# Encryption
def MutiURAG_Enc(Public_Key, Message,EG_Generator):  
    Support = random_support_gen(w_2)
    R1 = homogenous_matrix_gen_with_support(n,N2,Support)
    E = homogenous_matrix_gen_with_support(N1,N2,Support)
    R2 = homogenous_matrix_gen_with_support(n,N2,Support)
   
    
    U = R1 + Public_Key[0].transpose() * R2
    
    VV = E + Public_Key[1].transpose() * R2 
    V = Fold(Encoding_AG(Message,EG_Generator), N1, N2) + VV 
    
    ct = [U, V]
    return ct

# Decryption
def MutiURAG_Dec(Private_Key, Ciphertext, EG_Generator, r):  
    Fold_Noisy_Word = Ciphertext[1] - Private_Key[1].transpose() * Ciphertext[0]
    Noisy_Word = UnFold(Fold_Noisy_Word)
    return Decoding_AG(Noisy_Word, EG_Generator,r)


# Muti-UR-AG 
#(q,m,n,N1,N2,k,t,w_1,w_2) = (2,67,30,10,13,3,67,7,8) # Muti-UR-AG-128
#(q,m,n,N1,N2,k,t,w_1,w_2) = (2,83,38,12,14,3,83,8,9) # Muti-UR-AG-192
(q,m,n,N1,N2,k,t,w_1,w_2) = (2,113,45,13,15,3,113,9,10) # Muti-UR-AG-256


Fqm = GF(q**m)
Frob = Fqm.frobenius_endomorphism()
S = OrePolynomialRing(Fqm, Frob, 'x')

%time Public_Key, Private_Key = MutiURAG_KGen(q,m,n,N1,N2,k,t,w_1,w_2)

Message = random_vector(Fqm, k);  g = random_small_vec_gen(N1*N2, m)
%time Ciphertext = MutiURAG_Enc(Public_Key, Message, g)

r = w_1*w_2
%time Message_test = MutiURAG_Dec(Private_Key, Ciphertext, g, r)

# check correctness
Message_test == Message
