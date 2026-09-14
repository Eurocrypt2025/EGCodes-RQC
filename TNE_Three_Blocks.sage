# The Blockwise RD Problem with Three Blocks
# Theoretical Non-zero Equations (TNE)
# Theoretical Maximal Equations (TME) 

from math import comb
from itertools import product, combinations, combinations_with_replacement

def find_product(R,Re):
    results = []
    for prod in product(range(R + 1), repeat=Re):
        if sum(prod) == R: 
            results.append(prod)
    return results

def Theoretical_Nonzero_Equations(q,m,n,k,r,n1,n2,n3,r1,r2,r3):
    n0 = 0; r0 = 0; z0 = -(k+1); z1 = n1-(k+1); z2 = n1+n2-(k+1); z3 = n1+n2+n3-(k+1)
    Sum = 0
    if 0 < k+1 <= r1:
       results = find_product(r0-z0,3)
       for prod in results:
           Sum = Sum + comb(z1, r1+z0+prod[0])*comb(n2, r2+prod[1])*comb(n3, r3+prod[2])
       return Sum
    if r1 < k+1 <= n1:
       results = find_product(r1,3)
       for prod in results:
           Sum = Sum + comb(z1, prod[0])*comb(n2, r2+prod[1])*comb(n3, r3+prod[2])
       return Sum
    if n1 < k+1 <= n1+r2:
       results = find_product(r1-z1,2)
       for prod in results:
           Sum = Sum + comb(z2, r2+z1+prod[0])*comb(n3, r3+prod[1])
       return Sum
    if n1+r2 < k+1 <= n1+n2:
       results = find_product(r1+r2,2)
       for prod in results:
           Sum = Sum + comb(z2, prod[0])*comb(n3, r3+prod[1])
       return Sum
    if n1+n2 < k+1 <= n1+n2+n3:
       return comb(z3, r)


#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,15,1,6,5,5,5,2,2,2) 
# TNE & TME:  1100 & 1716;  0 < k+1 <= r1
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,15,2,6,5,5,5,2,2,2) 
# TNE & TME:  700 & 924;  r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,15,3,6,5,5,5,2,2,2) 
# TNE & TME:  400 & 462;    r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,15,4,6,5,5,5,2,2,2) 
# TNE & TME:  200 & 210;   r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,15,5,6,5,5,5,2,2,2) 
# TNE & TME:  84 & 84;   n1 < k+1 <= n1 + r1
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,15,6,6,5,5,5,2,2,2) 
# TNE & TME:  28 & 28;    n1 < k+1 <= n1 + r1
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,15,7,6,5,5,5,2,2,2) 
# TNE & TME:  7 & 7; n1 + r1 < k+1 <= n1 + n2

# (q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,15,6,6,6,5,4,2,2,2) 
# TNE & TME: 28 & 28   n1 < k+1 <= n1 + r2
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,15,7,6,7,5,3,2,2,2) 
# TNE & TME: 7 & 7     n1 < k+1 <= n1 + r2
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,17,8,6,8,5,4,2,2,2) 
# TNE & TME: 28 & 28   n1 < k+1 <= n1 + r2
(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,17,8,6,8,6,3,2,2,2) 
# TNE & TME: 25 & 28   n1 < k+1 <= n1 + r2
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,17,8,6,8,7,2,2,2,2) 
# TNE & TME: 15 & 28   n1 < k+1 <= n1 + r2

#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,20,7,6,14,3,3,2,2,2)  
# TNE & TME: 172 & 924   r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,20,8,6,14,3,3,2,2,2)  
# TNE & TME: 121 & 462   r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,20,9,6,14,3,3,2,2,2)  
# TNE & TME: 79 & 210    r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,20,10,6,14,3,3,2,2,2)  
# TNE & TME: 46 & 84    r1 < k+1 <= n1

#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,20,7,6,3,3,14,2,2,2)  
# TNE & TME: 924 & 924   
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,20,8,6,3,3,14,2,2,2)  
# TNE & TME: 462 & 462   
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,20,9,6,3,3,14,2,2,2)  
# TNE & TME: 210 & 210    
#(q,m,n,k,r,n1,n2,n3,r1,r2,r3) = (2,9,20,10,6,3,3,14,2,2,2)  
# TNE & TME: 84 & 84    



TME = comb(n-(k+1), r)
TNE = Theoretical_Nonzero_Equations(q,m,n,k,r,n1,n2,n3,r1,r2,r3)


print( )
print("Theoretical Non-zero Equations (TNE): ", TNE)
print("Theoretical Maximal Equations (TME) : ", TME)
print("TNE & TME:", TNE, "&", TME)
