# The Blockwise RD problem with Two Blocks
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

def Theoretical_Nonzero_Equations(q,m,n,k,r,n1,n2,r1,r2):
    n0 = 0; r0 = 0; z0 = -(k+1); z1 = n1-(k+1); z2 = n1+n2-(k+1)
    Sum = 0
    if 0 < k+1 <= r1:
       results = find_product(r0-z0,2)
       for prod in results:
           Sum = Sum + comb(z1, r1+z0+prod[0])*comb(n2, r2+prod[1])
       return Sum
    if r1 < k+1 <= n1:
       results = find_product(r1,2)
       for prod in results:
           Sum = Sum + comb(z1, prod[0])*comb(n2, r2+prod[1])
       return Sum
    if n1 < k+1 <= n1+n2:
       return comb(z2, r)


#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,1,4,6,6,2,2) 
# TNE & TME: 185 & 210; 0 < k+1 <= r1
#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,2,4,6,6,2,2) 
# TNE & TME: 120 & 126; r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,3,4,6,6,2,2) 
# TNE & TME: 70 & 70; r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,4,4,6,6,2,2) 
# TNE & TME: 35 & 35; r1 < k+1 <= n1
#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,5,4,6,6,2,2) 
# TNE & TME: 15 & 15; r1 < k+1 <= n1

#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,3,4,8,4,2,2) 
# TNE & TME: 53 & 70;  0 < k+1 <= r1
#(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,4,4,8,4,2,2) 
# TNE & TME: 31 & 35;  r1 < k+1 <= n1
(q,m,n,k,r,n1,n2,r1,r2) = (2,9,12,5,4,8,4,2,2) 
# TNE & TME: 15 & 15;  r1 < k+1 <= n1

TME = comb(n-(k+1), r)
TNE = Theoretical_Nonzero_Equations(q,m,n,k,r,n1,n2,r1,r2)

print( )
print("Theoretical Non-zero Equations (TNE): ", TNE)
print("Theoretical Maximal Equations (TME) : ", TME)
print("TNE & TME:", TNE, "&", TME)
