# Embedding degree

A subgroup G of an elliptic curve E(F_q) is said to have embedding degree k if the subgroup order divides q^k - 1, but does not divide q^i - 1 for 0 < i < k.

To have efficient pairings k needs to be small enough for computations in the field F_(q^k) to be feasible, but also big enough to have a required level of security.

Let E be an elliptic curve over the finite field F_q. Let N be the number of F_q points on E and l be a prime divisor of N. In crypto applications, E is always chosen so that N is divisible by a large prime l. Thus, solving dlog in the cyclic group of l points is hard.

## Torsion points

Let E be an elliptic curve over field K. Let l be a positive integer.

```
E[l] = {P from E(K~); l*P = 0} where K~ is algebraic closure of K
```

E[l] contains points with coordinates in K~, not just in K.

Let's say characteristic of K is not 2. What are the 2-torsion points of E: y^2 = x^3 + a*x + b? Obviously 2-torsion points are those where tangent is vertical. Tangent at (x1, y1) is:

```
(3*x1^2 + a) * x - 2*y1 * (y-y1) = 0
```

We can see that E[2] are points where y = 0. So we need to solve x^3 + a*x + b = 0 which gives us three points (and we also have point at infinity which is of torsion 2 as well). Are these three points pairwise different? If we take nonsingular curve, then yes. Recall that curve is singular iff its discriminant is 0 (discriminant of polynomial x^3 + a*x + b). But discriminat tells us: if it is 0, then polynomial has a multiple root.

So E[2] = {0, (e1,0), (e2,0), (e3,0)}. Later we will show that E[2] is isomorphic to Z_2 + Z_2.

Let's have a look at E[3] for curves in fields with characteristic different from 2 and 3. These are points P for which 2*P = -P which means that 2*P and -P have the same coordinate x.

Recall that for P1 = (x1, y1), P2 = (x2, y2), the coordinates of P3 = P1 + P2 are:

```
(x3, y3) = (m^2 - x1 - x2, m * (x1 - x2) - y1) where m = (y2-y1)/(x2-x1)
```

For P1 = P2:

```
(x3, y3) = (m^2 - 2*x, (x - x3) * m - y) where m = (3*x^2 + a) / 2*y
```

The equations x = m^2 - 2*x, m = (3*x^2 + a) / 2*y, and y^2 = x^3 + a * x + b give us:

```
m^2 * 4*y^2 = (3*x^2 + a)^2
3*x * 4*y^2 = (3*x^2 + a)^2
12*x*(x^3 + a*x + b) = (3*x^2 + a)^2
3*x^4 + 6*a*x^2 + 12*b*x - a^2 = 0
```

The discriminant of this polynomial is -6912(4*a^3 + 27*b^2)^2 which is nonzero, so we have 4 different zeros, meaning we have 8 different points in E[3] and also point at infinity, so 9 points. As we will see later:

```
E[3] = Z_3 + Z_3
```

What about for general n?

### E[n] is isomorphic to Z_n x Z_n

It can be shown (by division polynomials which enable us to see the degree of rational functions which determine endomorphism [n]) that endomorphism [n] (multiplication by n) is separable and of degree n^2.

For separable endomorhpisms it holds (Proposition 2.21 in Washington [2]): 

```
#ker(end) = deg(end) where end is endomorphism
```

That means that [n] has n^2 elements in the kernel, thus E[n] has n^2 elements. 

By structure theorem for finite abelian groups we know:

```
E[n] = Z_n_1 x ... x Z_n_k where n_i divides n_(i+1)
```

Let's take a prime l that divides n1. Then l divides all n_i. It holds: E[l] ⊆ E[n]. Because l|n_i, the order of E[l] is l^k. But we also know that order of E[l] is l^2, thus k = 2 and:

```
E[n] = Z_n_1 x Z_n_2
```

We know n^2 = n_1 x n_2, but we also know that n annihilates Z_n_1 x Z_n_2, thus n1|n and n2|n, so n_1 = n_2 = n and:

```
E[n] = Z_n x Z_n
```

## Balasubramanian-Koblitz theorem

Let E(F_q) contains a subgroup of order l where l does not divide q-1. Then for any positive integer k, E(F_q^k) contains all l^2 points of order l if and only if l divides q^k - 1.

Proof: It is known that if E(F_q^k) contains E[l] (all l^2 points of order l), then l divides q^k - 1 (where does it come from?). The conditions l divides #E(F_q) and l does not divide q-1 are not needed.

Conversely, suppose k > 1 and l divides q^k - 1.


page 48, pbc thesis
B-K paper


## Menezes-Okamoto-Vanstone (MOV) attack

MOV attack uses Weil pairing to move dlog problem in E(F_q) into F_(q^k). Dlog in finite fields can be attacked by index calculus methods and can be thus solved faster than dlog in elliptic curves.

We are trying to find x from x*P. If we have linearly independent P and Q, we have e(x*P, Q) = e(P, Q)^x where e(P, Q) is not 1. Thus, we can find x by breaking dlog in F_(q^k).

# Types of pairings (see [1])

Pairing is a function e:

```
e: G1 x G2 -> GT
```

where G1, G2, GT are cyclic groups of prime order l. If G1 = G2, we say it is a symmetric pairing.

Pairings can be separated in four different types:

 * G1 = G2
 * G1 != G2, but there is an efficiently computable homomorphism phi: G2 -> G1 (but not from G1 to G2)
 * G1 != G2 and there are no efficiently computable homomorphisms between G1 and G2
 * G2 is a non-cyclic group of order l^2.

Note that the cyclic groups of the same order are isomorphic, so there is an isomorhpism between G1 and G2, but it might not be efficiently computable. 

Let's say q is the size of the field over which elliptic curve E is defined. G1 is a subgroup of E(F_q), G2 is usually a subgroup of E(F_q^k), GT is a subgroup of F_q^k\*. Thus, there are three main parameters: field size q, embedding degree k, group size l.

## Type 1

When G1 = G2, Weil pairing becomes trivial since G1 is a cylic group and two elements are always dependent, let's say we have P, Q, then Q = l*P for some l. Then:

```
e(P, Q) = e(P, l*P) = e(P, P)^l = 1
```

In type 1 pairings we actually use a distortion map: 

```
psi: E(K)[r] -> E[r]\G1.
psi(x,y) = (-x, i*y) where i^2 = -1
```

## Type 2

Let's have a curve E over F_q with embedding degree k > 1. G1 is subgroup of E(F_q) of order l. For G2 we choose a random point from E(F_q^k)[l] and define `G2 = <Q>`. Q is published as system parameter.

Let's observe two maps - Frobenius map and trace map. Frobenius map is: pi: x -> x^q. Note that Frobenius map is a generator of Galois group Gal(F_q^k, F_q).

Trace map is:

```
Tr: E(F_q^k) -> E(F_q)
Tr(x) = x + pi^1(x) + pi^2(x) + ... + pi^(k-1)(x) = x + x^q + x^q^2 + ... + x^q^(k-1)
```

If we calculate Tr(x)^q:

```
Tr(x)^q = x^q + x^q^2 + ... + x = Tr(x)
```

Tr(x)^q = Tr(x) means Tr(x) is in F_q (because x -> x^q is an automorphism F_q^k -> F_q^k that preserves exactly F_q).

Furthermore, as Tr: E(F_q^k) -> E(F_q) is group homomorphism, if Q is from E[l] where l | #E(F_q), Tr(Q) is also in E[l] because l*Q = 0 and 0 = Tr(l*Q) = l*Tr(Q) and l is prime.


So we have a homomorphism from G2 to G1. The advantage of type 2 pairings is that any curve can be used and there is still a homomorphism from G2 to G1. On the other hand G2 does not have any special structure and it is impossible to sample randomly from G2 except by computing multiples of the generator, thus we cannot securely hash to G2. Also, membership testing can be a serious overhead.

The existence of homomorphism from G2 to G1 is sometimes require for security proofs. There are many schemes for which the security proof does not work if the system is implemented using type 3 pairings instead of type 2 pairings.

## Type 3

Let's have a pairing friendly curve E over F_q of embedding degree > 1. For G1 we take the subgroup of E(F_q) of order l. We can say G1 = ker(pi - [1]) ∩ E[l].

For G2 we take E[l] ∩ ker(pi - [q]). Note that [m] means x -> m*x in elliptic curve group. Thus, G2 are points on E[l] for which it holds:

```
(x^q, y^q) = q*(x,y)
```







[1] Galbraith, Steven D., Kenneth G. Paterson, and Nigel P. Smart. "Pairings for cryptographers." Discrete Applied Mathematics 156.16 (2008): 3113-3121.

[2] [9] Washington, Lawrence C. Elliptic curves: number theory and cryptography. CRC press, 2008.
