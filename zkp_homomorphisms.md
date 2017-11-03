# Zero-knowledge proofs of knowledge for homomorphisms using sigma_psi protocols

Most of the notes taken from [1].

Proofs of knowledge of a preimage under a homomorphism are a key building block in many constructions in applied cryptography. Essentially all constructions in applied (public key) cryptography are based on group homomorphisms. These are mappings psi: G -> H, where the domain is the group (G, +) and the co-domain is the group (H, *), such that:

```
psi(x1 + x2) = psi(x1) * psi(x2
```

In a proof of knowledge for a homomorphism, a prover demonstrates knowledge of a preimage x of y under psi (i.e., y = psi(x)). A proof of knowledge for a homomorphism is zero-knowledge, if a verifier learns nothing about a preimage of y under psi. If a prover succeeds in getting the verifier to accept with a probability larger than some threshold probability (the knowledge error), then the verifier can be asserted that the prover knows a solution (preimage). Knowing a solution means that an algorithm (the knowledge extractor) exists that, given the prover as a black-box, computes the desired solution.

All efficient zero-knowledge proofs of knowledge for homomorphisms are instance of the same protocol, which we will denote here as sigma_psi protocol. Well known example is Schnorr protocol (see [proofs1.md](https://github.com/miha-stopar/crypto-notes/blob/master/proofs1.md).

However, sigma_psi protocol is not efficient for all homomorphisms. It is efficient for example for Schnorr protocol, however for some homomorphisms, like exponentiation homomorphisms in hidden order groups (psi(x_1,...,x_n) = h_1^x_1 * ... * h_n^x_n), it is not. Exponetiation homomorhisms play an important role in applied cryptography (anonymous credentials, signatures, group signatures, voting systems, e-cash, multi-party computations), thus efficient zero-knowledge proofs are very desirable. In [1], efficient zero-knowledge proofs of knowledge for exponentiation homomorphisms in hidden order groups are designed.

## Efficiency limitations of sigma_psi protocols

The efficiency of a proof of knowledge (using sigma_psi protocol) for a given homomorphism psi is detrmined by the size of the knowledge error that can be achieved (the smaller the knowledge error can be made for psi, the more efficient proofs of knowledge can be obtained using sigma_psi protocol). For any psi there is a minimal knowledge error that is known to be achievable - let's call it minimal standard (MS) knowledge error.

It is important to note that the MS knowledge error is the smallest knowledge error that is known to be achievable using currently available knowledge extractors. That means there could be knolwedge extractors for the sigma_psi protocol that achieve a smaller knowledge error than the MS knowledge error, and hence would allow to overcome the efficiency limitations of the sigma_psi protocol.

## Efficient proofs of knowledge for exponentiation homomorphisms in hidden order groups

Results of [1] suggest that sigma_psi protocol cannot be used to obtain efficient proofs of knowledge for exponentiation homomorphisms. Thus, other techniques have to be used. In [1] three such techniques were proposed.

The first technique applies to exponentiation homomorphisms psi_E: Z^l -> H in hidden order groups H, provided that the prover (but not the verifier) knows the order of H. The underlying protocol is sigma_psi protocol in a novel setting where the common input additionally to psi_E and y contains a pseudo-preimage (v, w) of y under psi_E. A pseudo-preimage of y (where y from H) under psi: G -> H is a pair (v, w) with v from Z and w from G, such that y^v = psi(w). In this setting proofs of knowledge are obtained for essentially any psi_E in hidden order groups. The resulting proofs of knowledge are very efficient (arbitrary small knowledge error in single execution of the sigma_psi protocol).

The second technique is based on a novel protocol which will be here denoted by sigma_psi_plus. It yields so called proofs of knowledge in the auxiliary string model, which was introduced by Damgard [2]. The key idea is that an auxiliary string with a prescribed probability distribution is available to the prover and the verifier. The sigma_psi_plus protocol yields zero-knowledge proofs of knowledge in the auxiliary string model for any exponentiation homomorphism (thus also in hidden order groups). The resulting proofs of knowledge are very efficient when |image(psi_E)| contains only large prime factors. However, the setup protocol is rather inefficient (but in many cases the cost of the setup protocol does not matter).

The third techniqu is using sigma_psi_plus-WS protocol, which is a variant of the sigma_psi_plus protocol that removes the setup protocol and the associated cost. The sigma_psi_plus-WS protocol is a zero-knolwedge proofs of knowledge in random oracle model.

# Sigma_psi protocols

![sigma_psi_protocol](https://raw.github.com/miha-stopar/crypto-notes/master/img/sigma_psi_protocol.png)

TODO: start with 4.1, page 63




[1] Bangerter, Endre. Efficient zero knowledge proofs of knowledge for homomorphisms. Diss. Ruhr University Bochum, 2005.

[2] Ivan Damgard. Efficient concurrent zero-knowledge in the auxiliary string model. In Bart Preneel, editor, Advances in Cryptology — EUROCRYPT 2000, volume 1807 of Lecture Notes in Computer Science, pages 431–444. Springer Verlag, 2000.

