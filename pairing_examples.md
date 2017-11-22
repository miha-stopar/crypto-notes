
# Example of Weil pairing

Let's check elliptic curve E: y^2 = x^3 + x over F_59. Weil pairing on this curve is discussed in [1].

```
sage: E = EllipticCurve(GF(59), [1,0])
sage: E
Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 59
sage: E.cardinality()
60
```

We can see that it has 60 points.

Point at infinity (identity of the group):

```
sage: P = E(0); P
(0 : 1 : 0)
```

Let's have a look at r-torsion points. P is r-torsion point of E(K) when:

```
r * P = 0
```

Let's get 5-torsion points for E:

```
sage: P.division_points(5)
[(0 : 1 : 0), (25 : 29 : 1), (25 : 30 : 1), (35 : 28 : 1), (35 : 31 : 1)]
```

This points form a subgroup G of size 5. We denote it as E(F_59)[5].

It can be proved that for some integer k >= 1, E(F_(q^k))[r] contains exactly r^2 points and is the direct product of two cyclic groups of order r. For all k1 >= k it holds:

```
E(F_(q^1))[r] = E(F_(q^k))[r] 
```

Actually, r and q needs to be coprime, but this won't be discussed here.

We can thus denote E[r] = E(F_(q^k))[r].

We extend E now to E(F_(59^2)):

```
sage: E1 = E.base_extend(GF(59^2)); E1
Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 59^2
```

Let's observe 5-torsion points of E1:

```
sage: P1 = E1(0)
sage: P1.division_points(5)

[(0 : 1 : 0),
 (24 : 25*z2 + 17 : 1),
 (24 : 34*z2 + 42 : 1),
 (25 : 29 : 1),
 (25 : 30 : 1),
 (34 : 28*z2 + 45 : 1),
 (34 : 31*z2 + 14 : 1),
 (35 : 28 : 1),
 (35 : 31 : 1),
 (12*z2 + 53 : 5*z2 + 9 : 1),
 (12*z2 + 53 : 54*z2 + 50 : 1),
 (25*z2 + 9 : 16*z2 + 41 : 1),
 (25*z2 + 9 : 43*z2 + 18 : 1),
 (25*z2 + 17 : 7*z2 + 48 : 1),
 (25*z2 + 17 : 52*z2 + 11 : 1),
 (25*z2 + 25 : 29*z2 + 40 : 1),
 (25*z2 + 25 : 30*z2 + 19 : 1),
 (34*z2 + 34 : 16*z2 + 2 : 1),
 (34*z2 + 34 : 43*z2 + 57 : 1),
 (34*z2 + 42 : 7*z2 + 4 : 1),
 (34*z2 + 42 : 52*z2 + 55 : 1),
 (34*z2 + 50 : 29*z2 + 49 : 1),
 (34*z2 + 50 : 30*z2 + 10 : 1),
 (47*z2 + 6 : 5*z2 + 45 : 1),
 (47*z2 + 6 : 54*z2 + 14 : 1)]
```

There are 25 = 5^2, thus the whole E[5].

Weil pairing is function f: E[r] x E[r] -> F_(q^k).

In E(F_(59^2)) Weil pairing exists and is not trivial.

Because -1 is not quadratic residue in F_59 (there is no solution for x^2 % 59 = 58), we can extend F_59 with i: F_59[i]. It holds F_(59^2) = F_59[i] because two Galois fields of the same order are isomorphic.





[1] Lynn, Ben. On the implementation of pairing-based cryptosystems. Diss. Stanford University, 2007.


