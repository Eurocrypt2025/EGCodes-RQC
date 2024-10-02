# RQC-NIST

def random_small_vec_gen(n,t):
    B = matrix(Fqm.base_ring(),t,m,0)
    while B.rank() != t:
        B = matrix(Fqm.base_ring(),[vector(Fqm.random_element()) for i in range(t)])
    C = matrix(Fqm.base_ring(),n,m,0)
    while C.rank() != t:
        C = matrix(Fqm.base_ring(),n,t,[Fqm.base_ring().random_element() for _ in range(n*t)]) * B   
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
        C = matrix(Fqm.base_ring(),n,t,[Fqm.base_ring().random_element() for _ in range(n*t)]) 
    c = C * support
    return c

def random_small_vec_gen_with_support_one(n,support):
    t = len(list(support))
    support = support[0]**(-1) * support
    C = matrix(Fqm.base_ring(),n,t,0)
    while C.rank() != t:
        C = matrix(Fqm.base_ring(),n,t,[Fqm.base_ring().random_element() for _ in range(n*t)]) 
    c = C * support
    return c

def random_three_small_vec_gen_with_two_supports(n,support1,support2):
    c1 = random_small_vec_gen_with_support(n,support1)
    c3 = random_small_vec_gen_with_support(n,support1)
    t1 = len(list(support1)); t2 = len(list(support2))
    support = vector(list(support1) + list(support2))
    C2 = matrix(Fqm.base_ring(),n,t1+t2,0)
    while C2.rank() != t1+t2:
        C2 = matrix(Fqm.base_ring(),n,t1+t2,[Fqm.base_ring().random_element() for _ in range(n*(t1+t2))]) 
    c2 = C2 * support
    return c1,c2,c3

def rank_R(z):
    return matrix(Fqm.base_ring(),n,m,[vector(z[i]) for i in range(n)]).rank()

def support_R(z):
    return matrix(Fqm.base_ring(),n,m,[vector(z[i]) for i in range(n)]).row_space()

def support_list_R(z):
    C = matrix(Fqm.base_ring(),n,m,[vector(z[i]) for i in range(n)]).row_space().basis_matrix()
    t = C.nrows()
    return list(vector(Fqm,[C[i] for i in range(t)]))

def Encoding_Gabidulin(Message, Gabidulin_Support):
    f = S(Message.list())  # The message polynomial 
    return vector(f.multi_point_evaluation(Gabidulin_Support))

def Decoding_Gabidulin(Noisy_Word, Gabidulin_Support, r): 
    g_monomials = [Gabidulin_Support[i]**(q**j) for i in range(n) for j in range(k+r)] 
    SC2 = matrix(Fqm,n,k+r,g_monomials) 
    y_monomials = [Noisy_Word[i]**(q**j) for i in range(n) for j in range(r+1)] 
    SC1 = matrix(Fqm,n,r+1,y_monomials) 
    SC = block_matrix(Fqm,1,2,[SC1,SC2])
    #Solution = list(SC.right_kernel().random_element())  
    Solution = SC.right_kernel_matrix()[0].list()  
    V = S(Solution[0:r+1])
    N_vector = vector(Solution[r+1:k+2*r+1])
    N = S(list(N_vector))
    ff,re = N.left_quo_rem(-V)
    return vector(ff.list())
    
# Key Generation
def Blockwise_RQC_KGen(q,m,n,k,w,w1,w2):
    h = R.random_element()
    support = random_support_gen(w)
    x = R(list(random_small_vec_gen_with_support_one(n,support)))
    y = R(list(random_small_vec_gen_with_support_one(n,support)))
    s = x + h*y
    pk = [h,s]; sk = [x,y]
    return pk, sk

# Encryption
def Blockwise_RQC_Enc(Public_Key,Message,Gabidulin_Support): 
    support1 = random_support_gen(w1)
    support2 = random_support_gen(w2)
    vec_r1,vec_e,vec_r2 = random_three_small_vec_gen_with_two_supports(n,support1,support2)
    r1 = R(list(vec_r1)); r2 = R(list(vec_r2)); e = R(list(vec_e))
    u = r1 + Public_Key[0]*r2
    v = Encoding_Gabidulin(Message,Gabidulin_Support) +  vector(e + Public_Key[1]*r2)
    ct = [vector(u),v]
    return ct

# Decryption
def Blockwise_RQC_Dec(Private_Key,Ciphertext,Gabidulin_Support,r):  
    u = R(list(Ciphertext[0]))
    Noisy_Word = Ciphertext[1] - vector(u*Private_Key[1])
    return Decoding_Gabidulin(Noisy_Word, Gabidulin_Support,r)


# RQC-NIST
#(q,m,n,k,w,w1,w2) = (2,127,113,3,7,7,6) # RQC-NIST - 128
#(q,m,n,k,w,w1,w2) = (2,151,149,5,8,8,8) # RQC-NIST - 192
(q,m,n,k,w,w1,w2) = (2,181,179,3,9,9,7) # RQC-NIST - 256


Fqm = GF(q**m)
P1  = GF(q)['XX'].irreducible_element(n, algorithm = "minimal_weight")
P.<XX> = Fqm[]
R.<X> = P.quotient(P1)

Frob = Fqm.frobenius_endomorphism()
S = OrePolynomialRing(Fqm, Frob, 'x')

%time Public_Key, Private_Key = Blockwise_RQC_KGen(q,m,n,k,w,w1,w2)

Message = random_vector(Fqm,k); g = random_small_vec_gen(n,n)
%time Ciphertext = Blockwise_RQC_Enc(Public_Key,Message,g)

r = w*w1 + w2
%time Message_test = Blockwise_RQC_Dec(Private_Key,Ciphertext,g,r)

# check correctness
Message_test == Message
