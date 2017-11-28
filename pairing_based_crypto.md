# Pairing-based cryptography

Notes mostly taken from [1] and [2].

The following notation will be used:

```
k~: an algebraic closure of k
```

## Bilinear pairing or bilinear map

Let G_1, G_2 be two cyclic groups of the same prime order n. Let us view G_1 as an additive group with identity ∞ and G_2 as a multiplicative group with identity 1. Let P be a generator of G_1.

A mapping e: G_1 x G_1 -> G_2 is called bilinear map when:

 * e(aP, bQ) = e(P, Q)^(ab) for all P, Q from G_1 and a, b from Z_q* (bilinearity)
 * e(P, P) != 1, which means if P is a generator of G_1, e(P, P) is a generator of G_2 (non-degeneracy)
 * e can be efficiently computed (computability)

For all S, T from G_1 it holds:

 * e(S, ∞) = 1 and e(∞ ,S) = 1 
 * e(S, -T) = e(-S, T) = e(S, T)^(-1)
 * e(aS, bT) = e(S, T)^(ab) for all a, b from Z
 * e(S, T) = e(T, S)
 * if e(S, R) = 1 for all R from G_1, then S = ∞. 

A symmetric bilinear map is completely defined by the value it takes at e(P, P). There are n bilinear maps for G1, one of them is degenerate case. The other n-1 are equivalent up to a constant (e_1(S, T) = e_2(S, T)^c for all S, T from G1). What is difficult is to find a map that is efficiently computable.

## Discrete logarithm problem (DLP)

If we have an additively-written group G = <P> of order n and elements P and Q, DLP is the problem to find the integer x from [0, n-1] such that Q = xP. DLP is for example believed to be intractable in multiplicative group of a finite field and in group of points on elliptic curve over a finite field.

Due to the bilinearity property, DLP in G_1 can be efficiently reduced to the DLP in G_2.

This is because if we have (P, Q) in G_1 where Q = xP, then e(P, Q) = e(P, xP) = e(P, P)^x.

## Diffie-Hellman problem (DHP)

DHP is closely related to DLP. The problem is: given P, aP and bP, find abP.

The security of Diffie-Hellman key agreement protocol relies on DHP.

## Decisional Diffie-Hellman problem in cyclic groups with bilinear maps

Given P, aP, bP, cP, find if ab = c. As e(aP, bP) = e(P, P)^(ab), we can simply check whether this value is the same as e(P, P)^c. 

Thus, the cryptosystems relying on intractability of decisional Diffie-Hellman cannot be constructed in a cyclic group with bilinear pairing.

## Bilinear Diffie-Hellman problem (BDHP)

The security of many pairing-based protocols relies on the intractabilty of BDHP which is as follows: given P, aP, bP, cP, compute e(P, P)^(abc).

Hardnes of BDHP implies the hardness of DHP in both G_1 and G_2. The problem is generally assumed to be as hard as DHP in G_1 and G_2.

## Decisional Bilinear Diffie-Hellman Problem

Given P, aP, bP, cP, e(P, P)^w, determine whether w = abc.

## Some pairing-based protocols

### Three-party one-round key agreement

How can a key be exchanged between three parties in one round? Using Diffie-Hellman there are two rounds needed for three parties.

Using bilinear pairing, each of the three parties randomly selects a secret integer from [1, n-1] and broadcast the point to the other two parties (say Alice chooses a and broadcasts aP, Bob has b and broadcasts bP, Chris has c and broadcasts cP). Now each one of them can calculate a shared secret: Alice calculates it as e(bP, cP)^a, Bob as e(aP, cP)^b, Chris as e(aP, bP)^c. The shared secret is e(P, P)^(abc).

For parties that do not know either a, b, or c, the shared secret can be computed only as solving BDHP.

### Short signatures

Boneh, Lynn and Shacham (BLS) signature is created as follows:

 * Alice randomly chooses a from [1, n-1] which is her private key, A = aP is public key
 * H is hash function from {0, 1}* to G1\{∞}
 * signature on message m is created as S = aM where M = H(m)

The verifier just needs to check whether e(P, aM) = e(aP, M) (verifier knows S = aM, A = aP).

### Identity-based encryption

In identity-based cryptography a public key of Alice consists of her identifying information ID_A, like her e-mail address. A Trusted Third Party (TTP) use its private key to generate Alice'\'s private key from ID_A and securely transmit it to Alice.

Bob can encrypt messages for Alice using only ID_a and the TTP public key (Bob can encrypt message even before the private key for Alice is generated).

Such a system can be achived using scheme proposed by Boneh and Franklin:

 * we have bilinear pairing e for which BDHP is intractable
 * we have two hash functions H1: {0, 1}* -> G1\{∞} and H2: G_2 -> {0, 1}^l, where l is the bitlength of the plaintext
 * TTP private key is a randomly selected integer t from [1, n-1], its public key is T = tP
 * when Alice requests her private key d_A, TTP creates Alice''s identity string ID_A, computes d_A = tH_1(ID_A) and securely delivers d_A to Alice
 * Bob encrypts the message m from {0, 1}^l as Q_A = H1(ID_A), selects a random integer r from [1, n-1], computes R = rP, and finally computes c = m xor H2(e(Q_a, T)^r) and trnasmits the ciphertext (R, c) to Alice
 * Alice decrypts m = c xor H2(e(d_A, R)) (this works because e(d_A, R) = e(tQ_A, rP) = e(Q_A, tP)^r = e(Q_A, T)^r

Note that this scheme is secure against eavesdroppers, but is not resistant to chosen-ciphertext attacks. Let us say we have a ciphertext (R, c) and can somehow obtain a decryption m1 for (R, c1) where c1 is the same as c, only the first bit is flipped. Now we know m is the same as m1, only the first flip needs to be flipped. However, resistance to chosen ciphertext attacks can be achieved by using two additional hash functions which we will not cover here.

## Extension of pairing definition

Symmetric pairings can be instantiated by using suitable supersingular elliptic curves. To allow a wider range of curves to be used, the definition must be modified. Boneh, Lynn and Shacham used in [3] the following definition (multiplicative notation is used here).

Let G1, G2, G3 be cyclic groups of prime order r. Assume the DHP is hard in G1. Let phi: G2 -> G1 be an efficiently computable group isomorphism. Let g2 be a generator of G2 and g1 = phi(g2) (means g1 generates G2). A bilinear pairing e is an efficiently computable functions e: G1 X G2 -> G3 such that e(g1, g2) != 1 and e(g1^a, g2^b) = e(g1, g2)^(ab) for all a, b from Z.

This is sometimes called asymmetric pairing. It allows pairings to be constructed on ordinary curves. However, there is no known method to hash to an element of G2 such that its discrete log to some fixed base is unknown. We can still hash to elements of G1 and G3, and perform all other cryptographically useful operations in G2 (picking randome element, multiplications, inversions).

Due to a problem with hashing, the definition is further loosened.

### General bilinear pairing

Let r be a prime. Let G1, G3 be cyclic groups of order r. Let G2 be a group where each element has order dividing r (G2 is not necessarily cyclic). A bilinear pairing e is an efficiently computable function e: G1 X G2 -> G3 such that:

 * e(g1, g2) = 1 for all g2 from G2 iff g1 = 1, and similarly e(g1, g2) = 1 for all g1 from G1 iff g2 = 1 (nondegeneracy)
 * for all g1 from G1 and g2 from G2: e(g1^a, g2^b) = e(g1, g2)^(ab) for all a, b from Z (bilinearity).

Using this definition, the hardness assumptions need to be altered. For different schemes the assumptions (of problem difficulty) might be used in both G1 and G2, or a combination of the two. For example, we may need to assume that given g1, g1^x from G1 and g2 from G2, there is no efficient algorithm to compute g2^x (see [3]).

If r is not prime (composite group order is useful for some cryptosystems [4]), we must be aware that in this case all nondegenerate pairings are not equivalent up to a constant.

### Exponentiation as a bilinear pairing

Given some integer r > 1, let us take G1 = G3 = Z_r\* and G2 = Z_(r-1)\+. Define e: G1 X G2 -> G3 by e(g, a) = g^a.

According to the loosened definition above, this function is bilinear pairing.

## Divisors

Before we see how bilinear pairings are constructed, we first need to have a look at divisors. While rational functions of one variable are defined by zeros and poles, rational functions of two variables are defined by divisors.

Let's have an elliptic curve E: Y^2 = X^3 + A*X + B and a rational function of two variables f(X, Y). Then the numerator of f vanishes on some points of E (zeros) and denominator vanishes on some points of f (poles). Divisor of f is (it is only a notation):

```
div(f) = sum_(all points P of E) n_P * [P] where n_P is multiplicity at point P (poles have negative multiplicity)
```

More generally, we define a divisor on E to be any formal sum:

```
D = sum_(P from E) n_P * [P] with n_P from Z and n_P = 0 for all but finitely many P
```

We further define a degree and sum of a divisor:

```
deg(div(f)) = sum_(all points P of E) n_P
sum(div(f)) = sum_(all points P of E) n_P * P
```

Note that rational functions on E can have different expressions, for example f(x,y) = p(x,y)/q(x,y) = p1(x,y)/q1(x,y). See more about it in [algebraic_geometry.md](https://github.com/miha-stopar/crypto-notes/blob/master/algebraic_geometry.md))

Let's check an example E: y^2 = x^3 - x and the function f(x,y) = x/y on this curve. We see that:

```
x/y = y/(x^2-1)
```

Function f has thus two representations. We say that function is regular at point P if it has representation p/q where q(P) != 0. Thus f is regular at (0,0).


Example: y^2 = x^3 + 1 in GF(13)

```
sage: E = EllipticCurve(GF(13), [0, 1])
sage: E
Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 13
sage: E.points()
[(0 : 1 : 0), (0 : 1 : 1), (0 : 12 : 1), (2 : 3 : 1), (2 : 10 : 1), (4 : 0 : 1), (5 : 3 : 1), (5 : 10 : 1), (6 : 3 : 1), (6 : 10 : 1), (10 : 0 : 1), (12 : 0 : 1)]
```

Let's observe function g(x,y) = x^2 / y. So g(x,y) = p(x,y) / q(x,y) where p(x,y) = x^2 and q(x,y) = y. 

Let's homogenize E and g:

```
E: y^2 * z = x^3 + z^3
g(x,y,z) = x^2 / (y * z)
```

The zeros of p(x,y,z): (0 : 1 : 1), (0 : 12 : 1).
The zeros of q(x,y,z): (4 : 0 : 1), (10 : 0 : 1), (12 : 0 : 1), (0 : 1 : 0).

Note that (0 : 1 : 0) is point at infinity.

To know the divisor, we need to know the multiplicities (order) of these points. We need some [algebraic_geometry.md](https://github.com/miha-stopar/crypto-notes/blob/master/algebraic_geometry.md) here.

For a domain R that is not a field (because field doesn't have proper ideals) it is the same (see Fulton [6], section about Discrete Valuation Rings):

 * R is Noetherian and local, and the maximal ideal is principal.
 * There is an irreducible element t from R such that every nonzero z from R may be written uniquely in the form z = u * t^n, where u is a unit in R (invertible), t is a generator of maximal ideal, n is called order, t is called uniformizing parameter for R.

What is R in our case? It is O_P(V) - rational functions that are defined in P. The functions that are 0 at P do not have an inverse. V is a variety (in this case elliptic curve). While function field k(V) is obviously a field, O_P(V) is a ring.

It turns out O_P(V) is Noetherian (every ideal finitely generated), local (only one maximal ideal), and the maximal ideal is principal (generated by one element). That means every function z from O_P(V) can be written as z = u * t^n where u does not vanish at P, t vanish at P, and n is called order of the zero.

Let's have a look E: y^2 = x^3 + x at point P = (0, 0). The ideal (rational functions that vanish at P) of O_P(E) is generated by one function. Which? Obviously, functions x and y generate this ideal (all functions that vanish at (0, 0) are something multiplied by x or/and y). On E it holds x = y^2 * 1/(x^2 + 1), thus x is generated by y^2. Thus, y is a generator of the ideal.

So, what is divisor of g(x,y,z) = x^2 / (y * z) on E: y^2 * z = x^3 + z^3 in GF(13)? We saw that the zeros of g are (0 : 1 : 1), (0 : 12 : 1) and poles are (4 : 0 : 1), (10 : 0 : 1), (12 : 0 : 1), (0 : 1 : 0).

Let's check the function x^2 and the order of zero (0 : 1 : 1): it is obvious it is 2 (x is uniformizing parameter, unit u = 1). The same for (0 : 12 : 1).

What about poles (zeros of function y * z)? At (4 : 0 : 1) the uniformizing parameter is y, unit is z, order is 1. The same for (10 : 0 : 1) and (12 : 0 : 1). At (0 : 1 : 0) the uniformizing parameter is z, unit is y, order is 1.

```
div(g) = 2 * (0 : 1 : 0) + 2 * (0 : 12 : 1) - (4 : 0 : 1) - (10 : 0 : 1) - (12 : 0 : 1) - (0 : 1 : 0)
```

Theorem:

 * Let f and f1 be rational functions on E. If div(f) = div(f1), then there is a nonzero constant c such that f = c * f1.
 * Let D = sum_(P from E) n_P * [P] be a divisor on E. Then D is the divisor of a rational function on E iff: `deg(D) = 0 and sum(D) = 0`.

## Weil pairing

The set of points of order m (torsion points):

```
E[m] = {P from E: m * P = 0}
```

It holds (Silverman [8] Corollary III.6.4): Let E be an elliptic curve over F_p and assume that p does not divide m. Then there exists a value of k such that:

```
E(F_p^(j*k))[m] = Z_m x Z_m for all j >= 1
```

If m is prime and if K is a field, then we may view E[m] as a 2-dimensionl vector space over Z_m. Even if m is not prime, there is still {P1, P2} such that each P from E[m] can be written P = a * P1 + b * P2 where a, b are unique from Z_m.


Weil pairing takes as input a pair of points P, Q from E[m] and gives as output an m-th root of unity. Let's denote mi_m = {x from K~; x^m = 1}.

```
e_m : E[m] x E[m] -> mi_m
```

Weil pairing is a major tool in the study of elliptic curves - it can be for example used to prove Hasse's theorem on the number of points on an elliptic curve over a finite field.

Definition (see for example section 5.8.3 in [5]): Let be P, Q from E[m]. Let f_P and f_Q be rational functions on E satisfying:

```
div(f_P) = m * [P] - m * [0]
div(f_Q) = m * [Q] - m * [0]
```

The Weil pairing of P and Q is the quantity:

```
e_m(P, Q) = ( f_P(Q + S) / f_P(S) ) / ( f_Q(P - S) / f_Q(-S) )
```

where S from E is any point satisfying S not from {0, P, -Q, P-Q}.

Now, is Weil pairing well-defined - that means independent from the choice of f_P, f_Q and auxiliary point S? 

We can see that it is independent of the choice of f_P and f_Q because the rational functions on elliptic curve (see Theorem above) are by divisor defined up to a constant and in the Weil pairing the constants cancel each other out.

We can see that it is independent of the choice of S if have a look at the pairing as the function of S and compute its divisor:

```
F(S) = ( f_P(Q + S) / f_P(S) ) / ( f_Q(P - S) / f_Q(-S) )
```

It turns out the divisor is 0 and from the Theorem it follows, F(S) is constant.

Theorem (Washington [9], Theorem 11.2): Let E be an elliptic curve. Let D be a divisor on E with deg(E) = 0. Then there is a function f on E with div(f) = D iff sum(D) = 0.








[1] A. Menezes, “An Introduction to Pairing-Based Cryptography,” in 1991 Mathematics Subject Classification, Primary 94A60, 1991.

[2] Lynn, Ben. On the implementation of pairing-based cryptosystems. Diss. Stanford University, 2007.

[3] D. Boneh, B. Lynn, and H. Shacham. Short signatures from the Weil pairing. In Asiacrypt, volume 2248 of Lecture Notes in Computer Science, pages 514+, 2001.

[4] D. Boneh, E.-J. Goh, and K. Nissim. Evaluating 2-DNF formulas on ciphertexts. In Theory of Cryptography ’05, volume 3378 of Lecture Notes in Computer Science, pages 325–341. Springer-Verlag, 2005.

[5] Hoffstein, Jeffrey, et al. An introduction to mathematical cryptography. Vol. 1. New York: springer, 2008.

[6] Fulton, William. "Algebraic curves." An Introduction to Algebraic Geom (2008): 54.

[7] Shafarevich, Igor Rostislavovich, and Kurt Augustus Hirsch. Basic algebraic geometry. Vol. 2. Berlin: Springer-Verlag, 1994.

[8] Silverman, Joseph H. The arithmetic of elliptic curves. Vol. 106. Springer Science & Business Media, 2009.

[9] Washington, Lawrence C. Elliptic curves: number theory and cryptography. CRC press, 2008.

