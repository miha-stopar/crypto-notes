Notes mostly taken from [1].

# Decisional Diffie-Hellman (DDH) assumption

One of the most important cryptographic assumptions is DDH. It says: ([a], [r], [ar]) is computationally indistinguishable from ([a], [r], [z]) where a, r, z are uniform elements from Z_q and [a] means a*g (written additively) or g^a (written multiplicatively) and g is an element of a cyclic group of prime order q.

In bilinear groups DDH is no longer true - [ar] is no longer pseudorandom since e([a], [r]) = e(1, [ar]).

Boneh, Boyen, and Shacham [2] defined an alternative decisional assumption for G1 - given u, v, h, z from G1, (u, v, h, u^a, v^b, h^(a+b)) is computationally indistinguishable from (u, v, h, u^a, v^b, h^z). This is called (two-) linear decision assumption.

It can be written also: ([a1], [a2], [a1*r1], [a2*r2], [r1+r2]) is computationally indistinguishable from ([a1], [a2], [a1*r1], [a2*r2], [z]). It can be generalized to k-lin assumption.




[1] Escala, Alex, et al. "An algebraic framework for Diffie–Hellman assumptions." Journal of cryptology 30.1 (2017): 242-288.

[2] Boneh, Dan, Xavier Boyen, and Hovav Shacham. "Short group signatures." Annual International Cryptology Conference. Springer, Berlin, Heidelberg, 2004.
