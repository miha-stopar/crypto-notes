# sage -ipython pohlig-helllman.py

from sage.all import *
from random import randint
from multiprocessing import Process, Queue

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

def iterate(f, queue):
    a = randint(0, n-1)
    b = randint(0, n-1)
    point = None
    coordinates = (a,b)
    while True:
        point, coordinates = f(coordinates)
	x = point.xy()[0] 
        last_five_bits = ZZ(x) & 31
	# distinguishing property:
	if last_five_bits < 2:
    	    queue.put((point, coordinates))

if __name__ == '__main__':
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
    queue = Queue()

    processes = [Process(target=iterate, args=(f,queue,)) for x in range(4)]
    for p in processes:
        p.start()
    d = {}
    while True:
	point, coordinates = queue.get()
        if point in d:
	    coordinates1 = d[point]
	    if coordinates != coordinates1:
	        print coordinates
	        print coordinates1
    	        a = coordinates[0] - coordinates1[0]
	        b = coordinates1[1] - coordinates[1]
	        # it holds: a*P = b*Q
	        l = a * inverse_mod(b, n) % n
		print l
		break
	else:
	    d[point] = coordinates
	print "++"
	print point

    for p in processes:
        p.terminate()






