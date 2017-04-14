# Groups, rings, finite fields

Group is an abstraction of addition and subtraction, but the operation is not necessarily commutative (commutative groups are named Abelian groups).

Ring is a group with an additional operation (multiplication), where division is not necessarily available. Also, the addition in ring is commutative (example: integers).

Field is a ring where division is available (example: rational numbers).

## Group theory

Some definitions and facts used frequently in cryptography:

 * for integers x,y: GCD(x, y) is the greatest common divisor of x, y (for example GCD(12, 18) = 6)
 * for all integers x, y there exist integers a, b such that a\*x + b\*y = GCD(x, y)
 * if GCD(x, y) = 1 we say that x and y are relatively prime.
 * Z_n = {0, 1, 2, ... , n-1}; if we have x >= n, we take x % n (x modulo n)
 * x in Z_n has an inverse (there exists y for which x\*y = 1 (mod n)) if and only if GCD(x, N) = 1
 * Z_n* is by definition a set of invertible elements in Z_n
 * for prime p: Z_p\* = Z_p \ {0}

Some famous theorems from legendary mathematicians, on which a lot of cryptograhy is based:

 * Fermat (Fermat's little theorem): if p is prime, for all x from Z_p\* it holds: x^(p-1) = 1 in Z_p
 * Euler: Z_p\* is a cyclic group - there exists g from Z_p\* such that {1, g, g^2,...,g^(p-2)} = Z_p\* (g is called generator of Z_p\*)
 * Lagrange: the order of every subgroup H of G divides the order of G (thus: for all g from Z_p\* the order of g divides the order of Z_p\*, which is p-1) 
 * Euler (generalization of Fermat's little theorem, used for example in RSA): for all x from Z_n*: x^|Z_n\*| = 1 in Z_n
 * Cauchy: if G is a finite group and p is prime number dividing the order of G, then G contains an element of order p (there is a subgroup of order p in G)

Pierre de Fermat (1607-1665)             |  Leonhard Euler (1707-1783) | Joseph-Louis Lagrange (1736-1813) | Augustin-Louis Cauchy (1789-1857) (pictures from Wikipedia)
:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:
![screenshot](https://raw.github.com/miha-stopar/crypto-notes/master/img-wikipedia/fermat.jpg) | ![screenshot](https://raw.github.com/miha-stopar/crypto-notes/master/img-wikipedia/euler.jpg) | ![screenshot](https://raw.github.com/miha-stopar/crypto-notes/master/img-wikipedia/lagrange.jpg) | ![screenshot](https://raw.github.com/miha-stopar/crypto-notes/master/img-wikipedia/cauchy.jpg) 

It is quite surprising how much cryptography is based on these theorems. Some examples are given below.

### Finding a random prime

For example if we need a prime p of length 1024 bits, we can do:

 * step 1: choose a random integer p from [2^(1024), 2^(1025)-1]
 * step 2: if 2^(p-1) = 1 in Z_p return p, else go to step 1 

Usually a few hundred iterations are needed. The output is not guaranteed to be a prime, but it holds Pr[p not prime] < 2^(-60). The intuition: we choose a number x smaller than p and test x^(p-1), if this is not 1 in Z_p, then p is a composite (follows from Fermat). If it is 1, then p is with some probability a prime... Usually some further tests (like Miller-Rabin) are applied to check the primality of p.

### Finding a generator of a subgroup

The following fact that follows from theorems above can be used for example when trying to find the generator of a subgroup:

 * x is a generator of Z_p\* if and only if x^((p-1)/q) != 1 in Z_p for all q that divides p-1

This holds because if we say there exists x and q such that q|p-1 and x^((p-1)/q) = 1 in Z_p, then this means the order of x is (p-1)/q. This means that x is not the generator of Z_p\* (it does not have an order of p-1). So by proof of negation we proved that if x is a generator, then for all q that divides p-1 it holds x^((p-1)/q) != 1 in Z_p. For the other way around let us say that x is not a generator. Then q exists such that q < p-1 and x^q = 1. It also holds q | p-1 (Lagrange). If we define q_1 = (p-1)/q, then: x^((p-1)/q_1) = x^q = 1. So, there exists q (named q_1) for which it holds x^((p-1)/q) = 1. By proof of negation we proved that if for all q that divides p-1 holds x^((p-1)/q) != 1, then x is a generator.

The generator a of a subgroup of order q in Z_p\* can be thus found using the following procedure:

 * step 1: select an element y from Z_p\* and compute g = y^((p-1)/q
 * step 2: if g = 1 go to step 1, otherwise return g

### Extended Euclidean algorithm

For a given a and b the extended Euclidean algorithm finds x and y such that a * x + b * y = gcd(a, b). 

This algorithm can be used to calculate inverse in Z_p\*. If we write p * x + b * y = gcd(p, b), we have gcd(p, b) = 1 because p is prime. We are looking for the inverse of b. Now we have p * x + b * y = 1. But if we apply modulo p, we get b * y = 1. Thus, when we get y with the Extended Euclidean algorithm, we actually get the inverse of b modulo p.

Extended Euclidean algorithm works on the principle (assume a > b):
gcd(a, b) = gcd(b, a % b). Recursively we can thus translate gcd(a, b) into much smaller numbers. Like: gcd(76, 32) = gcd(32, 12) = gcd(12, 8) = gcd(8, 4) ... gcd(4, 0) (obviously the equation does not hold when b value reaches 0). Once b reaches 0 (like 0 = 8 % 4), the Eucliedan algorithm says that the gcd is the last value a (that is 4).

The x and y coefficents can be computed this way - when you pass a newly calculated a, b into gcd function (let us call this some_a and some_b), you pass also the coefficients with which these two values can be expressed by a and b.
For the first iteration, you simply say a = a * 1 + b * 0 and b = a * 0 + b * 1, the (x, y) coefficients are in this case (1,0) and (0,1). So we call:

gcd(a, b, (1,0), (0,1))

Or when we are deeper into recursion:

gcd(some_a, some_b, (some_a_x, some_a_y), (some_b_x, some_b_y))

It holds:

 * some_a = a * some_a_x + b * some_a_y
 * some_b = a * some_b_x + b * some_b_y

Function gcd calculates new_some_a, new_some_b and coefficients for new some_b (for new_some_a it already knows):

 * k = some_a / some_b (calculate as a whole number, ignore the remainder)
 * r = some_a % some_b

It holds:
r = some_a - k * some_b

Express some_a and some_b with a and b:
r = some_a_x * a + some_a_y * b - (some_b_x * a + some_b_y * b) * k

Calculate coefficients for r (which is a new some_b):
r = (some_a_x - some_b_x * k) * a + (some_a_y - some_b_y * k) * b

 * new_some_a = some_b
 * new_some_b = r

Finally, gcd calls itself with the following parameters:
gcd(new_some_a, new_some_b, (some_b_x, some_b_y), (some_a_x-some_b_x*k1, some_a_y-some_b_y*k1))

At each point in recursion we can express the current a and current b with the original a and b.
Thus, once the current b reaches 0, we can express current a (the actual greatest common divisor) with a and b - just what is needed for RSA.

See the recursion steps below for a = 32 and b = 12:

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/extended_euclid.png

The greatest common divisor is 4 and can be expressed as 4 = 32 * (-1) + 12 * 3.

## Chinese Remainder Theorem (CRT)

Let n_1, ..., n_r be pairwise relatively prime numbers. Let a_1,...,a_r be arbitrary integers. Then there is an integer a such that:

```
a = a_i (mod n_i) for each i = 1,...,r
```

Furthermore, a is unique modulo n = n_1 * ... * n_r.

Value a can be computed as follows. We define `N_i := n/n_i`. N_i and n_i are coprime and thus we can compute using extended Euclidean's algorithm M_i and m_i such that `M_i * N_i + m_i * n_i = 1`. Then:

```
a = a_1 * M_1 * N_1 + ... + a_r * M_r * N_r
```

This means there is a one-to-one correspondence between elements from Z_n and tuples from Z_(n_1) x ... x Z_(n_r).
Even more, this correspondence preserves the additive and multiplicative structure - it is an isomorphism. Formally:

The correspondence f between Z_n and direct product of rings of integers modulo n_i Z_(n_1) x ... x Z_(n_r) is a ring isomorphism.

That means f is a homomorphism and bijective (the inverse is also homomorphism). Homomorphism means:

```
f(e1) = e2 (identity into identity)
f(a + b) = f(a) + f(b)
f(a * b) = f(a) * f(b)
```

Ring isomorphism means that both rings can be considered to be the same for all questions related to addition and multiplication.

That means if we want to compute something in Z_n, we can always compute it in all Z_(n_i) first and then we apply the isomorphism to get the result in Z_n.
Let's say we want to compute c(x) where x is from Z_n. If we are able to compute c_i = c(x) % Z_(n_i), it means we have:

```
c_i = c(x) % n_i for all i = 1,...,r
```

Which means:

```
c(x) = c_i (mod n_i) for all i = 1,...,r
```

And we can compute c(x) by CRT.

## Some basics about groups

### Congruences

Let's have the congruence a * x = b (mod m). Let d > 0 be GCD(a, m). Let's denote a1 = a/d and m1 = m/d.

The congruence `a * x = b (mod m)` has solutions if and only if (iff) d|b.

If d|b, then there are exactly d solutions. If x0 is a solution, then other solutions are x0 + m1, x0 + 2 * m1, ..., x0 + (d-1) * m1.

### Normal subgroup

Subgroup H of group G is normal when:

```
x * H = H * x for all x from G.
```

### Direct product

Let G1, G2 be groups. If we define operation on G1 x G2 as (a,b) * (c, d) = (a*c, b*d), we get a group. We call it an *external direct product* of G1 and G2.

The notion of direct product of groups is used to factor a group into a product of smaller groups. This way in some cases we can get a complete characterization of a certain type of group.

Let G1 and G2 be normal subgroups of group G. Then G is called the *internal direct product* of G1 and G2 iff G = G1 * G2 and G1 intersection G2 is identity. 

If G = G1 x G2 (external product), then G is internal product of (G1, 1) and (1, G2).

If G is internal direct product of G1 and G2, then G and G1 x G2 are isomorphic.

From CRT if we have Z_n where n = p1^n1 * p2^n2 * ... * pk^nk, then

```
Z_n ~= Z_p1^n1 x Z_p2^n2 x ... x Z_pk^nk (~= means isomorphic)
```
p1, p2, ..., pk are relatively prime, so:

```
phi(n) = phi(p1^n1) * phi(p2^n2) * ... * phi(pk^nk)
phi(n) = (p1^n1 - p1^(n1-1)) * (p2^n2 - p2^(n2-1)) * ... * (pk^nk - pk^(nk-1)) 
phi(n) = p1^n1 * (1 - 1/p1) *  p2^n2 * (1 - 1/p2) * ... * pk^nk * (1 - 1/pk)
phi(n) = n * (1 - 1/p1) * (1 - 1/p2) * ... * (1 - 1/pk)
```

### Power residue in direct product

If we have x^n = a (mod m), a is called n-th power residue.

Let's have G = G1 x G2 where G1 is of order m and G2 is of order n.

Then the m-th power residues of G are exactly G2 (and n-th power residues are G1). This is because all m-th power residues are of the form (1, G2) (for every x from G1: x^m = 1), which is isomorphic to G2.

### Cyclic groups

Cyclic group is usually denoted by <g>, where g is the generator. For example Z_7* = {1,3,3^2, 3^3, 3^4, 3^5} = {1,3,2,6,4,5}, where the generator is 3.

All subgroups of a cyclic group are cyclic. For each d which divides n (order of <g>) there exists exactly one subgroup of order d and it is <g^(n/d)>.

#### Product of two cyclic groups

The product of two cyclic groups is cyclic if and only if their orders are relatively prime.

### Fundamental theorem of finite abelian groups

The theorem says: every finitte abelian group is a direct product of cyclic groups of prime-power order. The number of terms in the product and the orders of the cyclic groups are uniquely determined by the group.

This means every finite abelian group G is isomorphic to a group of the form:

```
Z_p1^n1 x Z_p2^n2 x ... x Z_pk^nk
```

where the pi^ni are uniquely determined by G.

However, as the product of two cyclic groups of relatively prime orders is a cyclic group (see above), there are various ways to express a given finite abelian group as a direct product of cyclic groups. For example Z_180 (180 = 2^2 * 3*2 * 5):

```
Z_180 = Z_4 x Z_9 X Z_5
Z_180 = Z_4 x Z_3 X Z_3 x Z_5
Z_180 = Z_2 x Z_2 X Z_9 x Z_5
Z_180 = Z_2 x Z_2 X Z_3 x Z_3 x Z_5
```

## Quadratic residue

We say a from Z_n is quadratic residue (special type of power residue, see above) if there exists some x such that x^2 = a (mod n). Otherwise a is called quadratic nonresidue.

Given N = p*q where p and q are two primes, it is difficult to find whether a is quadratic residue or nonresidue modulo N. 
Quadratic residues mod p is cyclic (subgroup of a cyclic group).
Note that QR_N where N is the product of two primes is also cyclic (see Cyclic groups above).

Note that distinguishing a quadratic residue from a nonresidue modulo p where p is prime is not a problem (the same holds for finding a root). According to the Euler's criterion, a is a quadratic residue modulo p if and only if a^((p-1)/2) = 1 (mod p). As p is an odd prime (in cryptography large primes are used) and thus (p-1)/2 is an integer, we can calculate a^((p-1)/2) to see whether a is a quadratic residue or not.

Note that polynomial x^2 - 1 = 0 has no roots besides 1 and -1, thus a^((p-1)/2) can be either 1 or -1 (we know that a^(p-1) = 1 (mod p)). Thus, a is quadratic residue if and only if a^((p-1)/2) = 1 (mod p) and nonquadratic residue if and only if a^((p-1)/2) = -1 (mod p).

### Euler's criterion

If p is odd prime and a an integer, then there exists x such that x^2 = a (mod p) if and only if a^((p-1)/2) = 1 (mod p).

Let's check both directions. 

If there exists x such that x^2 = a (mod p), then a^((p-1)/2) = x^(p-1). Due to Fermat's little theorem, x^(p-1) = 1 (mod p). 

We have a^((p-1)/2) = 1 (mod p) and want to find x such that x^2 = a (mod p). Let's say g is the generator of Z_p. There exists r such that g^r = a (mod p). We know g^(r*(p-1)/2) = 1 (mod p), which (because g^(p-1) = 1 (mod p)) means (p-1)|r*(p-1)/2 and we can thus see that r is even.

Thus there exists integer k such that r = 2*k and `a = g^(2*k) = (g^k)^2 (mod p)`.

### Quadratic residues modulo p where p is prime

Exactly half of Z_p\* are quadratic residues. Why: x^2 = (-x)^2 (mod p), which means if we take the first (p-1)/2 elements of Z_p\*, we have all quadratic residues (and it can be quickly shown that all these elements are different):

```
QR_p = 1^2, 2^2, ..., ((p-1)/2)^2
```

### Quadratic residues when N is product of two safe primes

Let p = 2 * p1 + 1, q = 2 * q1 +1, and N = p * q where p, p1, q, q1 are all primes (p and q are safe primes).

How can we get a generator g of QR_N?
The order of Z_N* is 4 * p1 * q1, the order of QR_N is 2 * p1 * q1. So we are searching for an element of order 2 * p1 * q1.

It turns out we can quickly get a generator: we choose a random a from Z_n* and check that a != 0 (mod p), a != 1 (mod p), a != -1 (mod p), and a != 0 (mod q), a != 1 (mod q), a != -1 (mod q). 
If the conditions are not met we try another random a. When we get a for which the conditions hold, we define generator g = a^2 (mod n).

Why? We need to know that if a = 0 (mod p), then a = 0 (mod n) (this follows from CRT or in other words: ring isomorphism maps zero to zero), so a^2 = 0 (mod n) and this is not a generator.  

If a = 1 (mod p) or a = -1 (mod p), then a^2 = 1 (mod p). So, a^2 = 1 (mod n) (again CRT). This is also not a generator.

So the conditions above rule out the wrong candidates. But are all the other elements generators?

The elements in QR_N have one of the following orders: 2, p1, q1, 2 * p1, 2 * q1.
If a is of order 2 (a^2 = 1 (mod n)), then a = 1 (mod n) or a = -1 (mod n), but we ruled out these elements (actually we ruled out elements for which a = 1 (mod p) ..., but due to ring isomorphis ...).
If a is of order p1 (a^p1 = 1 (mod n)), then a^p1 = 1 (mod q) and q-1 divides p1 (little Fermat). However, q-1 = 2 * q1 and we know that q1 and p1 are primes, so q-1 cannot divide p1. 
We do similarly for other orders and thus show that the only elements that are not ruled out, are of order 2 * p1 * q1.

## Schnorr group

Discrete logarithm problem means: given h such that h = g^x (mod p), find x.

C. P. Schnorr presented Schnorr groups in a paper describing efficient signatures for smart cards [1]. The goal was to reduce the amount of computation the smart card (or any other device) has to perform to generate a signature. Since then Schnorr groups are used in numerous protocols (see [proofs1.md](https://github.com/miha-stopar/crypto-notes/blob/master/proofs1.md) and [proofs2.md](https://github.com/miha-stopar/crypto-notes/blob/master/proofs2.md)).

Now, discrete logarithm is not of the same difficulty in all groups. Let's consider Z_n*. If |Z_n\*| can be factored to "small" factors for which discrete logarithm can be solved, then Pohlig-Hellman algorithm can be used to solve the discrete logarithm in Z_n*.

Thus, we want group of a prime order. How can we generate it? Taking Z_p* where p is prime is not ok, since the order of this group is p-1.

However, if we have two primes p and q where p = r * q + 1 for some integer r, we can find an element of order q in Z_p*. We just take g = a^r % p for some a. If it is not 1, then g is of order q.

The subgroup H of Z_p* generated by g is of order q.

```
H = {1, g, g^2 % p, g^3 % p, ..., g^(q-1) % p}
```

Note that this way speedup is obtained for computing g^x because x is smaller than q. On the other hand security is not sacrificed because the restriction to have an order of g much smaller than p provides no advantage in any of the known algorithms to break a discrete algorithm.

However, q needs to be large enough to resist birthday attack.

Also, index-calculus algorithms cannot be applied directly in H - first the index-calculus algorithm would need to be applied in Z_p* to then compute logarithms in H.

Often 160-bit q is used with 1024-bit p, but the recommendation is to use 256-bit q with 2048-bit p.

### Discrete logarithm in Schnorr group

Just as an interesting fact:

We have two primes p and q such that p = q * r + 1. Subgroup H of Z_p\* of order q is called Schnorr group.

For each element h from H there is an element a from Z_p\* such that h = a^r (mod p) (see Cyclic group section above - the order of Z_p\* is p-1, the order of H is q, thus H = <g^((p-1)/q)> = <g^r>).

Let's choose u, v from Z_p\*. Let's define s = u^r (mod p) and t = v^r (mod p).

Elements s and t are from H. If we somehow solve DL for s and t (that means we find x sucht that s^x = t (mod p)), we solved DL also for u and v.

So we have x such that:

```
s^x = t (mod p)
u^(r * x) = v^r (mod p)
```

If there exists y such that u^y = v (mod p):

```
u^y = v (mod p)
u^(r * y) = v^r (mod p)
```

So:

```
u^(r * x) = u^(r * y) (mod p)
```

This means that r * y = r * x + k * (p-1) for some k (because u^(p-1) = 1 (mod p) due to Little Fermat theorem). 

So we have (note that r|(p-1):

```
y = x + k * (p-1)/r
```

But we don't know k.

## Finite fields

A field is a structure with two operations: addition and multiplication (not necessarily the addition and multiplication in a traditional sense), represented by + and *. For operation +, all the elements must form a commutative group (identity element is denoted by 0). For operation \*, all elements except 0 (identitiy for +) must also form a commutative group (identity element is denoted by 1). Also, the distributive identity must hold: a * (b + c) = a * b + a * c. Fields are abstractions of familiar number systems, such as rational number, real number, and complex numbers.

The **order** of a finite field is the number of elements in the field. Finite fields are also called Galois fields.

Examples of *infinite* fields are rational numbers, real numbers and complex numbers. However, cryptography makes use of *finite* fields. Évariste Galois who was born in 1811 and killed at the age of 20 in a duel (picture from Wikipedia) found out that for any prime integer p and any integer n, a unique field with p^n elements exists, usually denoted GF(p^n). That means that if we take two Galois fields with the same number of elements, say A and B, there is always a mapping function f: A->B for which it holds f(a+b) = f(a) + f(b) and f(a * b) = f(a) * f(b). Such mapping function is called isomorphism.

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img-wikipedia/galois.jpg

### Prime Galois fields

If p is a prime number, then integers modulo p (these are {0,1,2,...,p-1}) with addition and multiplication performed modulo p, is a finite field of order p.

Prime Galois fields are often denoted also as Z/pZ or Z_p\*. A lot of crypto is actually done in prime Galois fields.

It needs to be pointed out: p needs to be prime, otherwise not all numbers have a multiplicative inverse. The integers modulo m (where m is not prime), do not form a field, but only a ring.

### Galois fields of prime power

In GF(p^n) elements are not numbers, but equivalence classes of polynomials whose coefficients belong to GF(p).

How come GF(p^n) is then a field? For this to achieve the addition and multiplication need to be properly defined. We have to see GF(p^n) as a collection of p^n n-dimensional vectors. Each coordinate in a vector is from GF(p).

Let us see an example in GF(2^8) (fields of order 2^m are called binary fields or characteristic-two finite fields). The two elements 83 and 249 can be written as: 

```
83 = 2^6 + 2^4 + 2^1 + 2^0
249 = 2^7 + 2^6 + 2^5 + 2^4 + 2^3 + 2^0
```

Or also:

```
83 = 01010011
249 = 11111001
```

Addition is defined as addition modulo p (p = 2 in this example), applied on each coordinate (when p = 2 this is the same as XOR operation): 01010011 + 11111001 = 10101010

Multiplication is defined (or is best explained) using polynomials. The two numbers above can be represented as:

```
83 = x^6 + x^4 + x^1 + x^0
249 = x^7 + x^6 + x^5 + x^4 + x^3 + x^0
```

We then multiply the two polynomials and get:

```
(x^6 + x^4 + x^1 + x^0) * (x^7 + x^6 + x^5 + x^4 + x^3 + x^0) = x^13 + x^12 + 2*x^11 + 2*x^10 + 2*x^9 + 2*x^8 + 3*x^7 + 3*x^6 + 2*x^5 + 3*x^4 + x^3 + x^1 + x^0
```

Now we need to reduce the coefficients modulo p and we get (the coefficients need to be in GF(2)):

```
x^13 + x^12 + x^7 + x^6 + x^4 + x^3 + x^1 + x^0
```

Now we need to reduce this polynomial once again because the exponents are not in GF(2^8) - the highest coefficient should be at most 7. This is done by dividing with an irreducible polynomial (the polynomial that cannot be factored) and keeping the remainder.

There is not only one polynomial that can be chosen for reducing - the multiplication will depend on it, but all the fields constructed using the properly chosen reducing polynomial will be isomorphic. For GF(2^8) we can choose for example x^8 + x^4 + x^3 + x^1 + x^0 (this is what AES uses, see ciphers.md).

Let us sage (SageMath) for a little demonstration. Let us create GF(2^8):

```
sage: k.<a> = GF(2^8)
sage: k
Finite Field in a of size 2^8
```

Now let us iterate over the elements in this field (only the 20 out of 256 are pasted here):

```
sage: for i,a in enumerate(k): print i, a
0 0
1 a
2 a^2
3 a^3
4 a^4
5 a^5
6 a^6
7 a^7
8 a^4 + a^3 + a^2 + 1
9 a^5 + a^4 + a^3 + a
10 a^6 + a^5 + a^4 + a^2
11 a^7 + a^6 + a^5 + a^3
12 a^7 + a^6 + a^3 + a^2 + 1
13 a^7 + a^2 + a + 1
14 a^4 + a + 1
15 a^5 + a^2 + a
16 a^6 + a^3 + a^2
17 a^7 + a^4 + a^3
18 a^5 + a^3 + a^2 + 1
19 a^6 + a^4 + a^3 + a
```

The default reducing polynomial is (not the same as in AES):

```
sage: k.polynomial()
a^8 + a^4 + a^3 + a^2 + 1
```

We just need to type into sage the polynomial we want to reduce and it will do the work for us:

```
sage: a^13 + a^12 + a^7 + a^6 + a^4 + a^3 + a^1 + 1
a^7 + a^4 + 1
```

So, the final result is a^7 + a^4 + 1. The arithmetic behind is:

```
a^13 + a^12 + a^7 + a^6 + a^4 + a^3 + a^1 + 1 = a^13 + a^12 + a^6 + a^3 + a^1 + (a^7 + a^4 + 1)
```

Note the part that is not in the brackets (what is in the bracket is the actual result), this should be expressible as a multiplication of reducing polynomial and some other polynomial. And it is:

```
a^13 + a^12 + a^6 + a^3 + a^1 = (a^5 + a^4 - a - 2) * (a^8 + a^4 + a^3 + a^2 + 1) - 2*a^7 + 2*a^4 + 4*a^3 + 2*a^2 + 2*a + 2
```

In GF(2^8) we discard all even coefficients. Thus, a polynomial (a^13 + a^12 + a^7 + a^6 + a^4 + a^3 + a^1 + 1) is really a reducing polynomial times some other polynomial (which is a^5 + a^4 - a - 2) plus (a^7 + a^4 + 1). 
Obviously, we could type into sage directly:

```
sage: (a^6 + a^4 + a^1 + a^0) * (a^7 + a^6 + a^5 + a^4 + a^3 + a^0)
a^7 + a^4 + 1
```

Some more properties of GF(256) - as xor is used for addition, + and - is the same operation (the inverse operation of xor is xor). A few polynomials which are all the same in GF(256):

```
a^13-1-a
a^7 + a^2

sage: a^13-1+a
a^7 + a^2

sage: a^13+1+a
a^7 + a^2

sage: a^13+1-a
a^7 + a^2

```

To simplify the implementation of the multiplication, the following equation can be used:

```
x^8 mod m(x) = m(x) - x^8
```

In AES we have m(x) = x^8 + x^4 + x^3 + x + 1, so:

```
x^8 mod (x^8 + x^4 + x^3 + x + 1) = (x^8 + x^4 + x^3 + x + 1) - x^8 = x^4 + x^3 + x + 1
```

This is because (for the second equality see above why + and - are the same operation; the last equality holds because the polynomial with degree more than 7 in GF(256) is calculated modulo irreducible polynomial, which in AES is (x^8 + x^4 + x^3 + x^2 + 1), so we simply take the remainder there):

```
x^8 = (x^8 + x^4 + x^3 + x + 1) * 1 - (x^4 + x^3 + x + 1) = (x^8 + x^4 + x^3 + x + 1) * 1 + (x^4 + x^3 + x + 1) = x^4 + x^3 + x + 1

```

If we have a polynomial b_7 * x^7 + b_6 * x^6 + ... + b_0 and we want to multiply it by 2 (means multiplying by x), we get b_7 * x^8 + b_6 * x^7 + ... + b_x * x. If b_7 = 0, then this is a polynomial with degree less than 8 and this is a final result. If b_7 = 1, we have (apply the equation above on x^8):

```
x^8 + (b_6 * x^7 + ... + b_0 * x) = (x^4 + x^3 + x + 1) + (b_6 * x^7 + ... + b_0 * x) = (b_6 * x^7 + ... + b_0 * x) + (x^4 + x^3 + x + 1)
```

So, if we are multiplying a number b_7,b_6,...,b_0 by 2 (irreducible polynomial is 00011011):

```
b_7,b_6,...,b_0 * 00000010 = b6,b_5,...,b_0,0  if b_7 = 0
b_7,b_6,...,b_0 * 00000010 = b6,b_5,...,b_0,0 + 00011011  if b_7 = 1

```

All GF(256) multiplication operations can be transformed into multiplications by 2 (and one addition operation), thus we can repeatedly apply the process above and get the result, like:

```
x * 9 = ((x * 2) * 2) * 2 + x
```

If for example we want to calculate 169 * 2 in GF(256):

```
10101001 * 00000010 = 01010010 + 00011011 = 1001001 = 64 + 8 + 1 = 73
```

In the first equality we used the formula above for b_7 = 1. So 169 * 2 = 73 in GF(256). It is not uncommon to have a pre-prepared multiplication tables in the actual implementations.

## Subfields and extension fields

A subset K1 of a field K is a **subfield** of K if K1 is also a field with respect to the operations of K. We say K is an **extension** of K1. A field F_(p^m) has precisely one subfield of order p^l for each divisor l of m (the elements of this subfield are the elements a from F_(p^m) for which it holds: a^(p^l) = a).

## Field characteristic

**Characteristic** of K is the smallest number ch for which ch * 1 = 0 where 1 is multiplicative identity element and 0 is additive identity (how many times we need to 1 + 1 + ... + 1 to get 0). If the sum is never 0, we say that the field has characteristic 0. 

Rational numbers, real numbers, complex numbers all have characteristic 0. If p is prime, the finite field GF(p^n) has characteristic p.

If K1 is a subfield of K, then K1 and K have the same characteristic.

### How AES makes use of Galois fields

Advanced Encryption Standard (AES) is a symmetric block cipher which makes use of GF(256) for some of its operations (SubBytes, MixColumns). For some explanation how AES works see ciphers.md.

# Problems without efficient solutions

Mathematical problems without known efficient solutions are the basis of cryptography.

## Factoring

Given N = p*q, where p and q are two primes of a similar size, find p and q (used for example in RSA).

## Discrete logarithm

Given h such that h = g^x (mod p), find x.

### Random self-reducibility for discrete logarithm

#### If we can find discrete logarithm for one base, we can find it for all bases

If we have an algorithm A which solves the discrete logarithm for g1 (that is we can always find x from g1^x (mod p)), we can find x from g^x (mod p) for any g.

This is because for a = g^x (mod p) we can find x1 such that a = g1^x1 (mod p). We can also find x2 such that g = g1^x2 (mod p). Now x1/x2 is the value we are searching for because:

```
g^(x1/x2) = (g1^x2)^(x1/x2) = g1^x1 = g (mod p)
```

#### If we can find an algorithm to break discrete logarithm for some fraction of inputs, we can break it for all inputs

Let's say we have a cyclic group G = <g>. If we have an algorithm A that solves discrete logarithm in G for a fraction 1/P of all inputs, then we actually have an algorithm to solve all inputs.

Let's say we want to find the discrete log of h. We try for various random values b to find a discrete log of h1 = h \* g^b by algorithm A. When we find b for which A was successful, we have log h = log h1 - b.

This is because:

```
g^(log h1 - b) = g^(log h1) / g^b = (h * g^b) / g^b = h
```

## RSA problem

RSA algorithm (see public_key_security.md) raises a message to an exponent, modulo a composite number n whose factors are not known:

```
c = m^e % n
```

Public key is a pair (n, e), private key is a pair (n, d) where d is inverse of e modulo phi_n (the number of invertible elements in Z_n). Message m is then calculated as (this holds due to Euler's theorem, see public_key_security.md):

```
m = c^d % n
```

RSA problem is to find m without knowing the private key which is the same as finding the e-th root of some arbitrary number c modulo n where n is is a composite number of unknown factorization (if n is prime we know how to compute e-th root of c).

## Flexible RSA problem (strong RSA assumption)

Strong RSA assumption differs from the RSA assumption in that the adversary can select the public exponent e. Given modulus n and a ciphertext c, adversary search for any plaintext m and public exponent e >= 3 such that c = m^e %n.

This is easier than solving the RSA problem, but still no efficient algorithm is known to solve it. Strong RSA assumption is the basis for many different cryptographic constructions.

Note that if we find v, x, y such that gcd(x, y) = 1 and the equation below holds, means that we solved flexible RSA problem:

```
v^x = u^y (mod n)
```

This is because we can find a and b such that:

```
a*x + b*y = gcd(x, y) = 1
```

And then:

```
v^x = u^y
v^((1-by)/a) = u^y
v^(1-by) = u^(ay)
v = u^(ay) * v^(by)
v^(1/y) = u^a * v^b
```

And:

```
m := v^(1/y)
m^x = u
```

## Distinguishing quadratic residue from quadratic nonresidue

Given N = p*q where p and q are two primes, find whether a is quadratic residue or nonresidue modulo N. 

Note that distinguishing a quadratic residue from a nonresidue modulo p where p is prime is not a problem (the same holds for finding a root). 

## Computational Diffie-Hellman (CDH) problem

Given some group G and group elements g, g^a, g^b, compute the value g^(ab).

Note that if we are able to find a discrete logarithm, we can immediately solve CDH.

## Decisional Diffie-Hellman (DDH) problem

Given some group G and group elements g, g^a, g^b, g^c, determine whether g^c = g^(ab).

If we are able to solve CDH, we can quickly solve DDH. If the other way around is possible, is not known.

An equivalent formulation of DDH is: given some group G and group elements g1, g2, g1^r1, g2^r2, determine whether r1 = r2.

Why is this equivalent? 

Let's name the first definition as DDH1 and the second as DDH2.

DDH1: given g, g^a, g^b, g^c, determine whether g^c = g^(ab).
DDH2: given g1, g2, g1^r1, g2^r2, determine whether r1 = r2.

### DDH1 => DDH2

We are able to solve DDH1 and we have g1, g2, g1^r1, g2^r2.
We would like to see whether r1 = r2.
It suffices to see whether g2^r1 = g2^r2.

There exists x such that g2 = g1^x (and we don't need to know it).
So g2^r2 = g1^(x * r2).
Our g^a is g1^r1, our g^b is g1^x = g2, our g^c is g1^(x * r2) = g2^r2.
Because we can solve DDH1, we can tell whether g2^r2 = g1^(r1 * x), so we can tell whether g2^r2 = g2^r1.

### DDH2 => DDH1

We are able to solve DDH2 and we have g1, g1^a, g1^b, g1^c.
We would like to see whether g1^(a * b) = g1^c.
So we need to see whether a * b = c.

Let's translate to DDH2. Let's define g2 as g1 = g2^a.
Now we have g1, g2, g1^c, g2^(a * b) = g1^b.
Because we can solve DDH2, we can tell whether c = a * b.

## Decisional Composite Residuosity Assumption (DCRA)

Let n = p * q be a product of two large primes. DCRA states that given n and integer z, it is hard to decide whether z is n-residue modulo n^2, meaning whether there exists y from Z_n^2* such that:

```
z = y^n (mod n^2)
```

The assumption was introduced by Paillier [2] in 1999.

The set of n-th residues is a multiplicative subgroup of Z_n^2 of order phi(n).

The problem of deciding n-th residuosity is denoted by CR[n]. Like the problem of deciding quadratic or higher degree residuosity, CR[n] is a random-self-reducible problem - all of its instances are polynomially equivalent. Each case is thus either uniformly intractable or uniformly polynomial.

## Composite Residuosity Class Problem

Like DCRA, Composite Residuosity Class Problem was presented in Paillier [2]. 

For w from Z_n^2*, find x from Z_n such that:

```
w = g^x * y^n mod n^2
```

where y is from Z_n*.


[1] C. P. Schnorr. Efficient Identification and Signatures for Smart Cards. In Crypto ’89, LNCS 435, pages 235–251. Springer-Verlag, Berlin, 1990.
[2] P. Paillier, Public-key cryptosystems based on composite residuosity classes, Advances in Cryptology — EUROCRYPT ’99, LNCS, vol. 1592, Springer Verlag, 1999, pp. 223–239.

