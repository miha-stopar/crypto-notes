# Embedding degree

A subgroup G of an elliptic curve E(F_q) is said to have embedding degree k if the subgroup order divides q^k - 1, but does not divide q^i - 1 for 0 < i < k.

To have efficient pairings k needs to be small enough for computations in the field F_(p^k) to be feasible, but also big enough to have a required level of security.

Let E be an elliptic curve over the finite field F_q. Let N be the number of F_q points on E and l be a prime divisor of N. In crypto applications, E is always chosen so that N is divisible by a large prime l. Thus, solving dlog in the cyclic group of l points is hard.

## Menezes-Okamoto-Vanstone (MOV) attack

MOV attack uses Weil pairing to move dlog problem in E(F_q) into F_(q^k). Dlog in finite fields can be attacked by index calculus methods and can be thus solved faster than dlog in elliptic curves.

We are trying to find x from x*P. If we have linearly independent P and Q, we have e(x*P, Q) = e(P, Q)^x where e(P, Q) is not 1. Thus, we can find x by breaking dlog in F_(q^k).

# Type of pairings

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







