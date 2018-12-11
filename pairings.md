## Briefly about endomorphisms

Quick notes on endomorphisms (Washington [2], section 2.9 or MIT lectures [7]).

Definition:
Morphism of projective curves is a rational map that is defined everywhere.

Theorem:
If C1 is a smooth projective curve then every rational map from C1 to a projective curve C2 is a morphism (proof in [8], see II.2.1).

Theorem:
A morphism of projective curves is either surjective or constant.

Definition:
An isogeny phi: E1 -> E2 of elliptic curves defined over k is a surjective morphism of curves that induces a group homomorphism E1(k~) -> E2(k~).

This definition is actually stronger than necessary. It could be: An isogeny phi: E1 -> E2 of elliptic curves defined over k is a non-constant rational map that sends point at infinity to point at infinity.

This is because any morphism of abelian varieties that preserves the identity element induces a group homomorphism.

Definition:
A morphism from an elliptic curve E/k to itself that fixes the distinguished point is called an endomorphism. An endomorphism that is also an isomorphism is an automorphism.

Except for the zero morphism, every endomorphism is an isogeny. The endomorphisms of an elliptic curve have a natural ring structure.


It holds:

 * each endomorphism can be written as (r1(x), r2(x) * y) where r1, r2 rational functions
 * def: endomorphism is separable if the derivative of r1(x) is not identically zero 
 * def: degree of endomorphism is max{deg(p(x), q(x)) where r1(x) = p(x)/q(x)
 * for each separable endomorphism end: #ker(end) = deg(end)




## Embedding degree

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

It turns out that this subgroup is precisely the kernel of trace map (see [13][12]). The kernel of the trace map can be obtained by a function P -> k*[P] - Tr(P). It can be easily seen that Tr(k*[P] - Tr(P)) = 0. This function is sometimes called anti-trace map.

# Pairings - examples in Sage

Mostly from [12] which is a great introduction to pairings.

Let's have q = 7691, E: y^2 = x^3 + 1, and F(q^2) constructed as F_q(u) where u^2 + 1 = 0 (Example 4.0.1 from [12]).

```
sage: q = 7691
sage: E = EllipticCurve(GF(q), [0, 1]); E
Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7691
sage: K.<x> = GF(q^2, modulus=x^2+1); K
Finite Field in x of size 7691^2
sage: E_ext = E.base_extend(K); E_ext
Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in x of size 7691^2
sage: E.cardinality()
7692 (note, this is 2^2 * 3 * 641)
sage: E_ext.cardinality()
59166864 (note, this is 2^4 * 3^2 * 641^2)
```

Now we will observe two points of order r = 641 and their Weil pairing:

```
sage: r = 641
sage: P = E_ext(2693, 4312)
sage: Q = E_ext(633*x + 6145, 7372*x+109)
sage: P.order()
641
sage: Q.order()
641
sage: P.weil_pairing(Q, 641)
6744*x + 5677
```

The result of Weil pairing lies in F(q^2) and is r-th root of unity: e(P, Q)^r = 1. Weil pairing is bilinear map e: G1 x G2 -> GT where all three groups are of order r. Bilinearity means: e(a*P, b*Q) = e(P, Q)^(a*b).

Note that Weil pairing is useful for cryptography only when r is huge, meaning that discrete logarithm cannot be computed in r-order subgroup.

Note also that the parameters above have been carefully chosen - r needs to divide #E(F(q^2)) so that E(F(q^2)) contains all r^2 r-torsion points; P and Q are r-order points (not in the same subgroup).

However, what is this mysterious Weil pairing? To construct Weil pairing we first need two rational functions of some special format. Rational function is any function p(x)/q(x) where p, q are polynomials. What about the special format? Actually, the divisor of a rational function needs to be of a special format. Divisor is something that tells us the zeros (zeros of polynomial p) and poles (zeros of polynomial q) of the rational function.

So, we are looking for two rational functions, let's say f_P and f_Q. f_P needs to have zero of degree r in point P and a pole of degree r in point at infinity, f_Q needs to have zero of degree r in point Q and a pole of degree r in point at infinity.

We can construct such functions by using Miller algorithm. Miller algorithm first finds a function that has zero of degree 2 in P, pole of degree 1 in P, and pole of degree 1 at infinty. This means the divisor is 2[P] - [P] - [0]. The algorithm then iteratively builds towards a function with divisor r[P] - [rP] - (r-1)[0] = r[P] - r[0] (note that rP = 0).

How do we get a function f such that div(f) = 2[P] - [P] - [0]?

Let's observe linear function (which is of course a rational function as well) e*x + f*y + g. This function has three zeros of degree 1 - it has three zeros because the zeros are the solutions of the following two equations:

```
e*x + f*y + g = 0
y^2 = x^3 + a*x + b
```

It can be shown that the discriminant of the polynomial which we get when we use the first equation in the second is nonzero which means that the three zeros are different.

## Torsion points

Example 4.1.1 from [12]:

```
sage: q=11
sage: K.<x> = GF(q^2, modulus=x^2+1);K
Finite Field in x of size 11^2
sage: E = EllipticCurve(GF(q), [0,4]);E
Elliptic Curve defined by y^2 = x^3 + 4 over Finite Field of size 11
sage: E.points()
[(0 : 1 : 0), (0 : 2 : 1), (0 : 9 : 1), (1 : 4 : 1), (1 : 7 : 1), (2 : 1 : 1), (2 : 10 : 1), (3 : 3 : 1), (3 : 8 : 1), (6 : 0 : 1), (10 : 5 : 1), (10 : 6 : 1)]
sage: len(E.points())
12

sage: E_ext = E.base_extend(K); E_ext
Elliptic Curve defined by y^2 = x^3 + 4 over Finite Field in x of size 11^2
```

Let's see for r = 3. We know there are r^2 (group of torsion points is isomorphic to Z_r x Z_r) torsion points, in our case 9. 

We know r does not divide q-1 and r divides q^2-1 = 120, so embedding degree k = 2. All torsion points can be found in E_ext (field extension from F(q) to F(q^2)). Let's see them:

```
sage: zero=E_ext(0)
sage: r = 3
sage: zero.division_points(r)
[(0 : 1 : 0),
 (0 : 2 : 1),
 (0 : 9 : 1),
 (8 : x : 1),
 (8 : 10*x : 1),
 (2*x + 7 : x : 1),
 (2*x + 7 : 10*x : 1),
 (9*x + 7 : x : 1),
 (9*x + 7 : 10*x : 1)]
```

There are 4 cyclic subgroups (there are always r+1 cyclic subgroups of r torsion points) of order 3 (at each subgroup there is also point at infinity which is omitted below):

```
{(0 : 2 : 1), (0 : 9 : 1)}
{(8 : x : 1), (8 : 10*x : 1)}
{(2*x + 7 : x : 1), (2*x + 7 : 10*x : 1)}
{(9*x + 7 : x : 1), (9*x + 7 : 10*x : 1)}
```

Let's observe Frobenius endomorphism:

```
sage: frob = K.frobenius_endomorphism()
sage: frob(9)
9
sage: frob(9*x+7)
2*x + 7
```

Frobenius endomorphism on torsion points:

```
(0 : 2 : 1) -> (0 : 2 : 1)
(0 : 9 : 1) -> (0 : 9 : 1)
(8 : x : 1) -> (8 : 10*x : 1)
(8 : 10*x : 1) -> (8 : x : 1)
(2*x + 7 : x : 1) -> (9*x + 7 : 10*x : 1)
(2*x + 7 : 10*x : 1) -> (9*x + 7 : x : 1)
(9*x + 7 : x : 1) -> (2*x + 7 : 10*x : 1)
(9*x + 7 : 10*x : 1) -> (2*x + 7 : x : 1)
```

If we take two different points from two different subgroups we will get two generators of all torsion points. Let's take p1 = (0 : 2 : 1) and p2 = (8 : x : 1). Each torsion point can be written as a * p1 + b * p2 where a, b from Z_r.

We can now write Frobenius endomorphism as 2x2 matrix with elements from Z_3.

We might for example take as generators of torsion group the following two elements:

```
e1 = (0 : 2 : 1)
e2 = (8 : x : 1)
```

or 

```
f1 = (2*x + 7 : x : 1)
f2 = (9*x + 7 : x : 1)
```

The matrix corresponding to e1, e2 is E = [[1, 0], [0, 2]]. The matrix corresponding to f1, f2 is E = [[0, 2], [2, 0]]. We can see that tr(E) = 1 + 2 = 3 = 0 and tr(F) = 0 + 0 = 0.




## Twisted curves

Twists of elliptic curves are used for efficient representations of elements in G2.

Let's say we have two curves in F_q:

```
E: y^2 = x^3 + a*x + b
E1: y^2 = x^3 + a*v^2*x + v^3*b
```

If we have v = c^2, we have a map psi E1(F_q) -> E(F_q):

```
psi(x, y) = (c^2 * x, c^3 * y)
```

Now, we want v to be a quadratic nonresidue so that we will be able to non-trivially extend the field F(q) to F(q^2).







Let's see an example:

```
sage: E = EllipticCurve(GF(19), [1, 6]); E
Elliptic Curve defined by y^2 = x^3 + x + 6 over Finite Field of size 19
sage: E1 = EllipticCurve(GF(19), [4, 10]); E1
sage: E.cardinality()
18
Elliptic Curve defined by y^2 = x^3 + 4*x + 10 over Finite Field of size 19
sage: E1.cardinality()
22

# now extend the field and observe curve in extended field F(19^2)
sage: E_ext = E.base_extend(GF(19^2)); E_ext
Elliptic Curve defined by y^2 = x^3 + x + 6 over Finite Field in z2 of size 19^2
sage: E1_ext = E1.base_extend(GF(19^2)); E1_ext
Elliptic Curve defined by y^2 = x^3 + 4*x + 10 over Finite Field in z2 of size 19^2
sage: E_ext.cardinality()
396
sage: E1_ext.cardinality()
396
```





# Barreto-Naehrig (BN) curves

BN curves [10] are curves of prime order (meaning that we take the whole curve, not only some prime order subgroup) and embedding degree k = 12. The equation of the curve is y^2 = x^3 + b with b != 0.

To recall - a subgroup G of an elliptic curve E(F_q) has an embedding degree or security multiplier k if the subgroup order r divides q^k - 1, but does not divide q^i - 1 for ) < i < k. Note that for pairings we need r|(q^k - 1) - pairing e: G1 x G2 -> GT where GT is subgroup of order r in F_q^k. 

The problem is to construct curves containing a subgroup with embedding degree k such that is it at once big enough to prevent the Frey-Ruck attack, but small enough that pairings are efficiently computable [9].

BN curves are constructed using some smart parameterisation. 

Let's have a look at BN256, here a prime p over which we form a basic field is of form 36*u^4 + 36*u^3 + 24*u^2 + 6*u + 1 where u = 1868033^3 = 6518589491078791937 (see for example Go implementation [11]). 

Order (number of elements in both, G1 and G2) is 36*u^4 + 36*u^3 + 18*u^2 + 6*u + 1 = 65000549695646603732796438742359905742570406053903786389881062969044166799969.

And t (trace) is 6*u^2 + 1.

Because t, n, p are parameterised, the space needed to store or transmit information about such curve is small.




[1] Galbraith, Steven D., Kenneth G. Paterson, and Nigel P. Smart. "Pairings for cryptographers." Discrete Applied Mathematics 156.16 (2008): 3113-3121.

[2] Washington, Lawrence C. Elliptic curves: number theory and cryptography. CRC press, 2008.

[3] Menezes, Alfred J., Tatsuaki Okamoto, and Scott A. Vanstone. "Reducing elliptic curve logarithms to logarithms in a finite field." iEEE Transactions on information Theory 39.5 (1993): 1639-1646.

[4] Schoof, René. "Nonsingular plane cubic curves over finite fields." Journal of combinatorial theory, Series A 46.2 (1987): 183-211.

[5] Schoof, René. "Elliptic curves over finite fields and the computation of square roots mod 𝑝." Mathematics of computation 44.170 (1985): 483-494.

[6] Miller, Victor. "Short programs for functions on curves." Unpublished manuscript 97.101-102 (1986): 44.

[7] https://math.mit.edu/classes/18.783/2017/LectureNotes5.pdf

[8] Joseph H. Silverman, The arithmetic of elliptic curves, Graduate Texts in Mathematics 106, second edition, Springer 2009.

[9] Barreto, Paulo SLM, Ben Lynn, and Michael Scott. "Constructing elliptic curves with prescribed embedding degrees." International Conference on Security in Communication Networks. Springer, Berlin, Heidelberg, 2002.

[10] Barreto, Paulo SLM, and Michael Naehrig. "Pairing-friendly elliptic curves of prime order." International Workshop on Selected Areas in Cryptography. Springer, Berlin, Heidelberg, 2005.

[11] https://github.com/golang/crypto/blob/master/bn256/constants.go

[12] http://www.craigcostello.com.au/pairings/PairingsForBeginners.pdf

[13] S. D. Galbraith. Pairings, volume 317 of London Mathematical
Society Lecture Notes, chapter IX, pages 183–213. Cambridge University
Press, 2005.
