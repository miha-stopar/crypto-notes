# sage -ipython pohlig-helllman.py

from sage.all import *
from random import randint

#E = EllipticCurve(GF(863056967), [334, 179])
E = EllipticCurve(GF(229), [1, 44])
#P = E.gens()[0]
P = E([5,116])
#Q = 4001*P
Q = E([155,166])

print "P: %s" % P
print "Q: %s" % Q

n = P.order()
print "P order: %s" % n

def H(X):
    x = X.xy()[0] 
    last_five_bits = ZZ(x) & 31
    return last_five_bits

def find_repetition(P, Q):
    a_list = []
    b_list = []
    for i in range(32):
        a_list.append(randint(0, n-1))
        b_list.append(randint(0, n-1))

    R = []
    for i in range(32):
        R.append(a_list[i] * P + b_list[i] * Q)

    def f((x1, x2)):
        X = x1 * P + x2 * Q
        j = H(X)
        a = a_list[j]
        b = b_list[j]
        return X + R[j], ((x1 + a_list[j]) % n, (x2 + b_list[j]) % n)

    count = 0
    tortoise, tortoise_coordinates, hare, hare_coordinates = P, (1,0), P, (1,0)
    while True:
        count += 1
        tortoise, tortoise_coordinates = f(tortoise_coordinates)
        hare, hare_coordinates = f(hare_coordinates)
        hare, hare_coordinates = f(hare_coordinates)
        if tortoise == hare:
	    if true:
		T = (1,0)
		for i in range(2*count):
    		    point, T = f(T) 
                    print point, T
    		    if i == count-1:
			print "----" 
  	    return tortoise_coordinates, hare_coordinates

def discrete_log(P, Q):
    while True:
        tortoise_coordinates, hare_coordinates = find_repetition(P, Q)
	if tortoise_coordinates != hare_coordinates:
    	    a = tortoise_coordinates[0] - hare_coordinates[0]
	    b = hare_coordinates[1] - tortoise_coordinates[1]
	    # it holds: a*P = b*Q
	    l = a * inverse_mod(b, n) % n
	    return l

l = discrete_log(P, Q)
print l





