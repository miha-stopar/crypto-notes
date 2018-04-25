# ZKPs

This tries to be a synopsis of ZKPs. Some more info (but also more scattered and not really organized) is in [proofs1.md](https://github.com/miha-stopar/crypto-notes/blob/master/proofs1.md).

There are two crucial properties about ZKPs:

 * proof of knowledge - there exists an algorithm which can extract the knowledge if the prover is used as a black-box and can be rewinded
 * zero knowledge - there exists a simulator which can simulate accepting transcripts which cannot be distinguished from the transcripts with a real verifier

In Schnorr protocol the extractor goes like (there is a secret s and publicy known h = g^s mod p):

 * prover outputs g^r mod p
 * simulator chooses challenge c1
 * prover returns z1 = r + c1 * s mod q
 * verifier checks whether g^z1 = g^r * (g^s)^c1 = g^r * h^c1
 * prover is rewinded so that it returns the same g^r and the protocol is repeated (the verifier checks whether g^z2 = g^r * (g^s)^c2 = g^r * h^c2)

Now the extractor divides both equations:

```
g^(z2-z1) = h^(c2-c1)
```

Note that this is happening in Schnorr group which is cyclic with order q:

```
g^(z2-z1) = h^(c2-c1) = g^(s*(c2-c1))
z2-z1 = s * (c2-c1) (mod q)
cdiff_inv * (z2-z1) = s (mod q)
```

Extractor now knows the secret: s = cdiff_inv * (z2-z1).

Whatever z1 and z2 the prover sent (potentially not of the form r + c*s), he obviously knows the value s for which g^s = h (mod p).

But note that the extractor needs to be able to compute cdiff_inv. For this the group order needs to be known (q in Schnorr group). Alternatively, the extractor works if (c1-c2) divides (z1-z2), which is the case in Damgard-Fujisaki commitments (where order is not known). Note also that if the challenge space is {0, 1}, the inverse is not needed (as c1 - c2 = 1), so in hidden order groups we have an extractor if we use binary challenges, however this is not efficient. The extractor shows that if the prover is able to succesffuly run the protocol two times, then it must know the secret (there is no way z1 and z2 can be prepared without knowing the secret).

In RSA group which is hidden order, the extractor still works due to the strong RSA assumption. We have g^(z2-z1) = h^(c2-c1) - due to the strong RSA assumption it holds that (z2-z1) = u * (c2-c1) for some u. So:

```
g^(u*(c2-c1)) = h^(c2-c1)
```

That doesn't necessarily mean that we have u such that g^u = h, it might be:

```
(g*b)^u = h 
```

where b^(c2-c1) = h. But it turns out b can only be 1 or -1.

In Schnorr protocol the simulator for proving the zero-knowledge property:

 * simulator chooses random c from challenge space and random z, it outputs transcript (g^z * (g^s)^(-c), c, z) which is accepting

However, these outputs cannot be distinguished from a real conversation between prover and verifier only if there challenges are chosen randomly as well - that means we trust verifier that it chooses challenges randomly. This kind of ZKP is called honest-verifier ZKP.

If we want to make Schnorr protocol to be ZKP (maybe it is anyway, but it is not known for sure), we can limit the challenges space to be small and execute protocol many times. In this case we can provide a simulator:

 * simulator chooses random c from challenge space and random z, it output the same transcript as above - if verifier accepts, the transcript is ok, otherwise simulator rewinds back the verifier and chooses other random c and provide the transcript - it repeats until the transcript is accepted (it won't take long as the challenge space is small). 
Because the knowledge error is 1/2 (if challenges space is {0, 1}), the protocol needs to be repeated multiple times.

An alternative to executing protocol many times with small challenges is approach proposed by Damgard [2] - here instead of g^r a commitment to it is sent and thus verifier cannot choose challenge based on the first message (challenge is thus chosen randomly).

# Zero-knowledge proofs of knowledge for homomorphisms using sigma_psi protocols

Most of the notes taken from [1]. The thesis proposes proofs of knowledge in hidden order groups (because Damgard-Fujisaki is only demonstration of a knowledge, not proof of knowledge).

Proofs of knowledge of a preimage under a homomorphism are a key building block in many constructions in applied cryptography. Essentially all constructions in applied (public key) cryptography are based on group homomorphisms. These are mappings psi: G -> H, where the domain is the group (G, +) and the co-domain is the group (H, *), such that:

```
psi(x1 + x2) = psi(x1) * psi(x2)
```

In a proof of knowledge for a homomorphism, a prover demonstrates knowledge of a preimage x of y under psi (i.e., y = psi(x)). A proof of knowledge for a homomorphism is zero-knowledge, if a verifier learns nothing about a preimage of y under psi. If a prover succeeds in getting the verifier to accept with a probability larger than some threshold probability (the knowledge error), then the verifier can be asserted that the prover knows a solution (preimage). Knowing a solution means that an algorithm (the knowledge extractor) exists that, given the prover as a black-box, computes the desired solution.

All efficient zero-knowledge proofs of knowledge for homomorphisms are instance of the same protocol, which we will denote here as sigma_psi protocol. Well known example is Schnorr protocol (see [proofs1.md](https://github.com/miha-stopar/crypto-notes/blob/master/proofs1.md).

However, sigma_psi protocol is not efficient for all homomorphisms. It is efficient for example for Schnorr protocol, however for some homomorphisms, like exponentiation homomorphisms in hidden order groups (psi(x_1,...,x_n) = h_1^x_1 * ... * h_n^x_n), it is not. This is because in a hidden order group it is hard to compute a non-zero multiple of a random group element.

Exponetiation homomorhisms play an important role in applied cryptography (anonymous credentials, signatures, group signatures, voting systems, e-cash, multi-party computations), thus efficient zero-knowledge proofs are very desirable. In [1], efficient zero-knowledge proofs of knowledge for exponentiation homomorphisms in hidden order groups are designed.

## Efficiency limitations of sigma_psi protocols

The efficiency of a proof of knowledge (using sigma_psi protocol) for a given homomorphism psi is detrmined by the size of the knowledge error that can be achieved (the smaller the knowledge error can be made for psi, the more efficient proofs of knowledge can be obtained using sigma_psi protocol). For any psi there is a minimal knowledge error that is known to be achievable - let's call it minimal standard (MS) knowledge error.

It is important to note that the MS knowledge error is the smallest knowledge error that is known to be achievable using currently available knowledge extractors. That means there could be knolwedge extractors for the sigma_psi protocol that achieve a smaller knowledge error than the MS knowledge error, and hence would allow to overcome the efficiency limitations of the sigma_psi protocol.

## Efficient proofs of knowledge for exponentiation homomorphisms in hidden order groups

Results of [1] suggest that sigma_psi protocol cannot be used to obtain efficient proofs of knowledge for exponentiation homomorphisms. Thus, other techniques have to be used. In [1] three such techniques were proposed.

The first technique applies to exponentiation homomorphisms psi_E: Z^l -> H in hidden order groups H, provided that the prover (but not the verifier) knows the order of H. The underlying protocol is sigma_psi protocol in a novel setting where the common input additionally to psi_E and y contains a pseudo-preimage (v, w) of y under psi_E. A pseudo-preimage of y (where y from H) under psi: G -> H is a pair (v, w) with v from Z and w from G, such that y^v = psi(w). In this setting proofs of knowledge are obtained for essentially any psi_E in hidden order groups. The resulting proofs of knowledge are very efficient (arbitrary small knowledge error in single execution of the sigma_psi protocol).

The second technique is based on a novel protocol which will be here denoted by sigma_psi_plus. It yields so called proofs of knowledge in the auxiliary string model, which was introduced by Damgard [2]. The key idea is that an auxiliary string with a prescribed probability distribution is available to the prover and the verifier. The sigma_psi_plus protocol yields zero-knowledge proofs of knowledge in the auxiliary string model for any exponentiation homomorphism (thus also in hidden order groups). The resulting proofs of knowledge are very efficient when |image(psi_E)| contains only large prime factors. However, the setup protocol is rather inefficient (but in many cases the cost of the setup protocol does not matter).

The third technique is using sigma_psi_plus-WS protocol, which is a variant of the sigma_psi_plus protocol that removes the setup protocol and the associated cost. The sigma_psi_plus-WS protocol is a zero-knolwedge proofs of knowledge in random oracle model.

# Sigma_psi protocols

![sigma_psi_protocol](https://raw.github.com/miha-stopar/crypto-notes/master/img/sigma_psi_protocol.png)

C = C(k) will denote a challenge set, k is security parameter.

Sigma_psi protocol for homomorphisms with a finite domain is described in a diagram above. It is generalization of Schnorr for arbitrary homomorphisms.

The sigma_psi protocol on a diagram is not defined for homomorphisms psi: G -> H with an infinite domain, such as exponetiation homomorphisms psi_E: Z^l -> H. The reason is that in the first step of the prover's computation the choice of a uniform random element from the infinite domain is not possible. However, for exponentiation homomorphismus we can restrict the domain of psi_E to a finite subset G of the domain Z^l. But one needs to be careful in order for the resulting protocol to be (honest-verifier) zero-knowledge. We need to introduce two finite subsets G, G1 of Z^l. Let's do only for l = 1.

```
delta_x = delta_x(k)
xc = xc(k)
```

Given a sequence of challenge sets C(k) we define:

```
gamma = max({abs(c): c from C(k)})
```

Using an auxiliary security parameter k_s = poly(k) we define:

```
G = {xc - delta_x, xc + delta_x} 
G1 = {-2^k_s * gamma * delta_x, 2^k_s * gamma * delta_x}
```

![sigma_psi_infinite_protocol](https://raw.github.com/miha-stopar/crypto-notes/master/img/sigma_psi_infinite_protocol.png)

Sigma_psi protocol for exponentiation homomorphisms with infinite domain is depicted above.

Theorem:

* sigma_psi_protocol for homomorphisms with a finite domain is honest-verifier zero-knowledge, and if for the challenge set C(k) we have |C(k)| <= poly(k), then it is zero-knowledge

* sigma_psi protocol for exponetiation homomorphisms is statistical honest-verifier zero-knowledge, and if for the challenge set C(k) we have |C(k)| <= poly(k), then it is statistical zero-knowledge

The first part is discussed in the description of Schnorr protocol in [proofs1.md](https://github.com/miha-stopar/crypto-notes/blob/master/proofs1.md).

For the second part we first prove the honest-verifier zero-knowledge property. We choose random c1 from C, random s1 from G1, set t1 = sigma_psi(s1)*sigma_psi(x)^(-c1)*sigma_psi(xc)^c1, and outputs (t1, c1, s1). We know that t1 is determined with s1, so we can only check that the transcripts (c, s) and (c1, s1) have statistically indistinguishable distributions. Because the verifier is assumed to be honest: c and c1 are equally distributed. We need to check s and s1.

We know that s is uniformly distributed in {-2^k_s * gamma * delta_x + c*(x-delta_x), 2^k_s * gamma * delta_x + c*(x-delta_x)} and s1 is distributed uniformly in G. The probability that s is outside G is less than 1/2^k_s. Thus the proof is honest-verifier statistical zero-knowledge.

The for perfect zero-knowledge is a combination of the above and the approach of proving that Schnorr is zero-knowledge for small challenge space.

## Proof of knowledge property

We discuss for which homomorphisms sigma_psi protocol is a proof of knowledge. For sigma_psi protocol we can extract pseudo-preimage of y (that means w such that psi(w) = y^v). This is pseudo-preimage extractability.

For knowledge extractor we first use pseudo-preimage extractor and then solve the pseudo-preimage (PP) problem: given a pseudo-preimage (v, w), compute x such that psi(x) = y.

Note that for homomorphism psi: G -> H and y from H, if (v, w) is a pseudo-preimage of y and gcd(v, |y|) = 1, then y = psi(v^(-1) * w).

Let's see why: 

```
y^v = psi(w)
gcd(v, |y|) = 1 => v * v^(-1) = 1 (mod |y|)
y^(v*v^(-1)) = y^(k*|y| + 1) = y
y = psi(w)^v^(-1) = psi(w * v^(-1))
```

Note that we need to know |y| to compute preimage.

When pseudo-preimage (PP) problem is solvable, we have a knowledge extractor for sigma_psi protocol.





TODO: 129




[1] Bangerter, Endre. Efficient zero knowledge proofs of knowledge for homomorphisms. Diss. Ruhr University Bochum, 2005.

[2] Ivan Damgard. Efficient concurrent zero-knowledge in the auxiliary string model. In Bart Preneel, editor, Advances in Cryptology — EUROCRYPT 2000, volume 1807 of Lecture Notes in Computer Science, pages 431–444. Springer Verlag, 2000.

[3] Ivan Damgard and Eiichiro Fujisaki. An integer commitment scheme based on groups with hidden order. In Advances in Cryptology — ASIACRYPT 2002, volume 2501 of LNCS. Springer, 2002.

