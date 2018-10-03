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

Recall: 

```
E[n] = E[F_q~][n] where F_q~ is closure of F_q
```

We will see in what follows when E(F_q) contains all n^2 n-torsion points of E(F_q~).

Lemma (see [3], Lemma 3): If gcd(n, q) = 1, then E[n] is contained in E(F_q) iff:

 * n^2 | #E(F_q)
 * n | q-1
 * TODO: check [4]

Let's have E(F_q) an elliptic curve over F_q, the finite field on q elements. Let q = p^m, where p is the characteristic of F_q.

By Hasse theorem, the order of E(F_q) is q + 1 - t where |t| <= 2*sqrt(q). E(F_q) is said to be supersingular if p divides t. It seems like Neal Koblitz was the first to conjecture that the supersingular curves have a small embedding degree.

From Hasse theorem we know t^2 <= 4*q. It actually turns out that t^2 can only be 0, q, 2*q, 3*q, or 4*q.

Moreover, the group structure of the supersingular curves is as follows:

 * if t^2 = q, 2*q, or 3*q, then E(F_q) is cyclic
 * if t^2 = 4*q, then either E(F_q) = Z_(sqrt(q)-1) ⊕ Z_(sqrt(q)-1) or E(F_q) = Z_(sqrt(q)+1) ⊕ Z_(sqrt(q)+1), depending on whether t = 2*sqrt(q) or t = -2*sqrt(q) respectively

The above are results from Schoof's paper [4] which answers on the following question: given a finite field F_q and integer N >= 0, how many projectively inequivalent non-singular elliptic curves are there over F_q that have exactly N points defined over F_q?

Weil pairing is a function:

```
e_n: E[n] x E[n] -> F_q~
```

Lemma (Lemma 5 in [3]): Let E(F_q) be an elliptic curve such that E[n] is contained in E(F_q), and where n is a positive integer coprime to q. Let P be from E[n] and of order n. Then for all P1, P2 from E[n], P1 and P2 are in the same coset of <P> within E[n] iff e_n(P, P1) = e_n(P, P2).

Proof:

If P1, P2 are from the same coset, then P1 = P2 + k*P and e_n(P,P1) = e_n(P,P2)*e_n(P,P^k = e_n(P,P2).

Conversely, suppose that P1, P2 are from a different coset. Let's choose Q such that (P,Q) generates E[n] = Z_n x Z_n. Then: P1 - P2 = a1*P + a2*Q where a2*Q is not identity.

```
e_n(P,P1) = e_n(P,P2 + a1*P + a2*Q) = e_n(P,P2) * e_n(P, P)^a1 * e_n(P, a2*Q) = e_n(P,P2) * e_n(P, a2*Q)
```

We are done if we prove e_n(P, a2*Q) is not identity. If e_n(P, Q) = 1, then for each element R = a*P + b*Q from E[n] we have e_n(P, a*P + b*Q) = e_n(P,P)^a * e_n(P,Q)^b = 1. But from non-deneneracy of e_n we then have that P is identity which is a contradiction.

## Structure of E(F_q)

Classical result from Max Deuring (1941):

E(F_q) is either cyclic or isomorphic to a product of two cyclic groups Z_n1 x Z_n2 where n2|n1.

Some other questions regarding the structure of E(F_q):

 * let e_q(E) be the size of the largest cyclic subgroup of E(F_q) (note that n1 = e_q(E)) - is it typical for E to have a large exponent e_q(E)?
 * how often is E cyclic?

### Number of points on E(F_q)

Number of points on E(F_q) can be computed in polynomial time using Schoof's algorithm [5].

Schoof's algorithm uses Frobenius endomorphism pi and its characteristic equation:

```
pi^2 - [t]*pi + [q] = 0
```

What does that mean? Let's use sage and curve y^2 = x^3 + x over finite field of size 59. The value t for this curve is 0, so characteristic equation is:

```
pi^2 + [q] = 0
```

Let's check:

```
E = EllipticCurve(GF(59), [1,0])
E1 = E.base_extend(GF(59^2))
q = 59
q2 = q^2

for p in E1.points():
    if p != E1(0):
        print(E1(p[0]^q2, p[1]^q2, 1) + q*p)
```

The output is (0 : 1 : 0) for all points.

Quick notes on endomorphisms (Washington [2], section 2.9):

 * each endomorphism can be written as (r1(x), r2(x) * y) where r1, r2 rational functions
 * def: endomorphism is separable if the derivative of r1(x) is not identically zero 
 * def: degree of endomorphism is max{deg(p(x), q(x)) where r1(x) = p(x)/q(x)
 * for each separable endomorphism end: #ker(end) = deg(end)

We can see that for pi, deg(r1) = 0, so Frobenius endomrphism is not separable. Another example of non-separable endomorphism is multiplication by p in characteristic p.

Now let's see how we get characteristic polynomial of Frobenius endomorphism. We can restrict each endomorphism to E[n] because endomorphism maps n-torsion points to n-torsion points. E[n] can be viewed as Z_n x Z_n and we can represent endomorphism with 2 x 2 matrix A. 

```
A = [[a b] [c d]]
```

We will need Cayley-Hamilton theorem. For a matrix A, we define its characteristic polynomial:

``` 
p(lambda) = det(A - lambda) where lambda is identity matrix multiplied by lambda
``` 

Cayley-Hamilton says p(A) = 0. In out case for 2 x 2 matrix A, this means:

```
p(lambda) = (a - lambda) * (d - lambda) - b*c = lambda^2 - lambda*(a + d) a*d - b*c = lambda^2 - lambda*tr(A) + a*d - b*c 
```

Note that trace of A (written as tr(A)) is the sum of diagonal elements.

By Cayley-Hamilton:

```
A^2 - A * tr(A) + det(A) = 0
```

By definition t = tr(A).




pi^n - 1 is separable (glej rjchen)
glej Washington 53



https://math.stackexchange.com/questions/1268144/question-about-characteristic-polynomial-of-the-frobenius-endomorphism-on-ellipt

Frobenius map is not separable:
https://people.cs.nctu.edu.tw/~rjchen/ECC2009/11_FrobeniusEndomorphism.pdf

http://pages.cs.wisc.edu/~cdx/MathJournal/October/CharPoly04.pdf


https://math.mit.edu/classes/18.783/2015/LectureNotes7.pdf




https://perso.univ-rennes1.fr/reynald.lercier/file/LM95.pdf
http://web.math.sinica.edu.tw/bulletin_ns/20144/2014406.pdf
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.370.3295&rep=rep1&type=pdf
https://math.mit.edu/classes/18.783/2017/LectureNotes9.pdf

https://crypto.stanford.edu/miller/miller.pdf

Also, n1 and n2 can be determined in probabilistic polynomial time by an algorithm due to Miller [6].


## Menezes-Okamoto-Vanstone (MOV) attack

MOV attack [3] reduces the elliptic curve logarithm in a curve E over a finite field F_q to the discrete logarithm in a suitable extension field F_q^k of F_q. This is done by establishing an isomorphism between <P> and the subgroup of n-th root of unity in F_q^k, where n denotes the order of P. The isomorphism is given by Weil pairing.

Dlog in finite fields can be attacked by index calculus methods and can be thus solved faster than dlog in elliptic curves.

For the MOV algorithm we need to be able to pick points P uniformly and randomly on E(F_q) in probabilistic polynomial time. This can be done by randomly choosing an element x1 from F_q. If x1 is x-coordinate of some point in E(F_q), then we can find y1 such that (x1,y1) from E(F_q) by solving a root finding problem in F_q. From Hasse theorem, the probability that x1 is the x-coordinate of some point in E(F_q) is at least (q+1-2*sqrt(q))/2*q, so at least 1/2 - 1/sqrt(q).



We are trying to find x from x*P. If we have linearly independent P and Q, we have e(x*P, Q) = e(P, Q)^x where e(P, Q) is not 1. Thus, we can find x by breaking dlog in F_(q^k).

## Balasubramanian-Koblitz theorem

Let E(F_q) contains a subgroup of order l where l does not divide q-1. Then for any positive integer k, E(F_q^k) contains all l^2 points of order l if and only if l divides q^k - 1.

Proof: It is known that if E(F_q^k) contains E[l] (all l^2 points of order l), then l divides q^k - 1 (where does it come from?). The conditions l divides #E(F_q) and l does not divide q-1 are not needed.

Conversely, suppose k > 1 and l divides q^k - 1.


page 48, pbc thesis
B-K paper



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

[3] Menezes, Alfred J., Tatsuaki Okamoto, and Scott A. Vanstone. "Reducing elliptic curve logarithms to logarithms in a finite field." iEEE Transactions on information Theory 39.5 (1993): 1639-1646.

[4] Schoof, René. "Nonsingular plane cubic curves over finite fields." Journal of combinatorial theory, Series A 46.2 (1987): 183-211.

[5] Schoof, René. "Elliptic curves over finite fields and the computation of square roots mod 𝑝." Mathematics of computation 44.170 (1985): 483-494.

[6] Miller, Victor. "Short programs for functions on curves." Unpublished manuscript 97.101-102 (1986): 44.
