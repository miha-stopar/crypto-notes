# execute: sage -ipython pohlig-helllman.py

from sage.all import *
import pdb

E = EllipticCurve(GF(863056967), [334, 179])
P = E.gens()[0]
Q = 4001*P

print "P: %s" % P
print "Q: %s" % Q

n = P.order()
print "P order: %s" % n
print "factorization of order of P: %s" % n.factor()

# l_1 = z_0 + z_1 * p_1 + ... + z_(e_1-1) * p_1^(e_1-1)

ls = []
for ind, factor in enumerate(n.factor()):
    # calculating l_ind
    print "----------------------------"
    p = factor[0]
    e = factor[1]
    print "p: %s, e: %s" % (p, e)
    zs = []
    for j in range(e):
	# calculating z_j
        smallQ = ZZ(n/(p**(j+1))) * Q
        smallP = ZZ(n/(p**(j+1))) * P
    	x = discrete_log(smallQ, smallP, operation='+')
	d = 0
	for index, i in enumerate(zs):
	    d += i * (p**index)
        z = (x - d)/(p**j) # z_j 
	zs.append(z)
    l_ind = 0
    for k in range(e):
	l_ind += zs[k] * p**k
    print "l_ with index %s: %s" % (ind, l_ind)
    ls.append(l_ind)
# Chinese Remainder Theorem:
fs = []
for i in n.factor():
    fs.append(i[0]**i[1])
print ls
print fs
l_reconstructed = crt(ls, fs)
print l_reconstructed





