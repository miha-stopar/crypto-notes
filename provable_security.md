# Provable security

The basis of any public-key cryptosystem is a one-way function - a function y = f(x) that is easy to evaluate but computationally infeasible to find the inverse x  f^(-1)(y). RSA for example uses y = x^e (mod n), where n is an integer whose prime factors are secret (we cannot compute inverse without knowing factors and factoring of large numbers is believed to be hard). 
Some public-key cryptosystems are based on function y = g^x where g is a fixed generator of a group of prime order in which discrete logarithm is believed to be hard to find.

Note: these notes are mostly taken while reading (somehow controversial, see Goldreich [2]) essay from Koblitz and Menezes [1].

## Rabin encryption

It was soon realised that breaking a public-key crypto system is not necessarily equivalent to solving the underlying mathematical problem. In 1979 designed an encryption function that could be proved to be invertible only if you know how to factor (large) n.

### Key generation

Choose two large distinct primes p and q. To simplify the computation of square roots modulo p and q use p = q = 3 (mod 4). Public key is n = p \* q, private key is (p, q).

When p = 3 (mod 4) the square root of c in mod p can be calculated

```
(c^((p+1)/4))^2 = c^((p+1)/2) = c^((p-1)/2) * c = c (mod p)
```

The last step follows from Euler's criterion (see [number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md)).

Hence the two square roots of c mod p are c^((p+1)/4) mod p and -c^((p+1)/4) mod p.

### Encryption

The message m is encrypted by `c = m^2 mod n`.
Note that four different values of m produce the same ciphertext.

### Decryption

Let's define a square root of c as root_c. We can compute root_c % p and root_c % q using the formula above (key generation). Let's say root_c_p is a square root of c in Z_p and root_c_q is a square root of c in Z_q. We have:

```
root_c_p = root_c % p
root_c_q = root_c % q
```

So:

```
root_c = root_c_p (mod p)
root_c = root_c_q % (mod q)
```

And we can compute root_c by CRT (see CRT section in number_theory.md). 
CRT says we first compute m_1 and m_2 such that `m_1 * p + m_2 * q = 1`. Then

```
root_c = root_c_p * m_2 * q + root_c_q * m_1 * p
```

However, we have two possible values for root_c_p (let's say root_c_p and -root_c_p) and root_c_q, which means we have four square roots of c.


### Provable security of Rabin encryption

If someone is able to compute two root_c values which are not additive inverses, for example:

```
root_c1 = root_c_p * m_2 * q + root_c_q * m_1 * p
root_c2 = root_c_p * m_2 * q - root_c_q * m_1 * p
```

then:

```
root_c1 + root_c2 = 2 * root_c_p * m_2 * q
```

We can see that q divides this value. What about p? We know that `m_1 * p + m_2 * q = 1. If p divides m_2, then we get `p * (m_1 + some_k * q) = 1`, but this is not possible as both factors are integers.
Also, p does not divide root_c_p, because that would mean that root_c_p is 0 (mod p). 
However, Rabin encryption does now work for messages that are 0 mod p or mod q (but these are very rare, because p and q are large). 

That means we can use gcd (which is fast) to compute q = gcd(root_c1 + root_c2, n). That means we are able to factor a large integer which we know it is not true.

But what if only one root is discovered? Now we need to remember the types of attack (see security_definitions.md).

Consider for example an adversary that is restricted to at most a known message attack (he knows signatures c_1, ..., c_r for some set of messages m_1, ..., m_r which was not chosen by him).
Let's say he knows the root of for example c_1. If he is able to find another of c_1, this second root is with probability 1/2 not an additive inverse of m_1.
Therefore, the adversary can with probability 1/2 factor n (the probability quickly gets higher when he tries different c_i). But in the known message attack, the adversary cannot choose messages to be decrypted, so he cannot get the second root of c_1 (if it is not in the initial set m_1, ..., m_r).

However (see [3]), if an adversary has the power of a chosen message attack, he can choose message m and asks for the signature of m^2 % n. With probability 1/2 he obtains a root that can help him to factorize n. Goldwasser, Micali and Rivest presented in [3] the signature scheme that is robust against an adaptive chosen message attack.

Note that plain RSA is also vulnerable to chosen-ciphertext attack - see section 'Chosen-ciphertext attack on plain RSA' above, although in this case only one message is decrypted, while in Rabin encryption private key is obtained and all subsequent messages can be decrypted.

One might think that chosen-ciphertext attack is too powerful to be worried about (who would give you a decryptions for arbitrary messages), but just think of Bleichenbacher attack for example [4].

## Some history of provable security

In 1980s it became clear (Rabin encryption ... ) that having a one-way function is not enough for security of public-key cryptography systems.

First, Goldwasser and Micali [5][6] introduced a notion of probabilistic encryption. Using a probabilistic encryption a given plaintext is never encrypted into the same ciphertext - examples are RSA PKCS1 and OAEP paddings.

Similarly, in Rabin encryption if you pad a message m with a random r before squaring mod n, the chosen-ciphertext attack is not possible anymore. But as mentioned in [3] also the proof of security does not hold anymore.

For probabilistic encryption Goldwasser and Micali [5][6] defined two notions of security for passive adversary (cannot get decryptions of chosen ciphertext) - semantic security (the adversary cannot obtain any information at all, this notion can be used for active adversary too) and indistinguishability (it has been extended for active adversary in [7], see definition at the top of this page).

In 2002 it has been proved [7][8] that indistinguishability and semantic security under chosen-ciphertext attacks are equivalent.
However, it turns out it is much easier to prove a public-key encryption scheme to be secure using indistinguishability than semantic security.
 
In 1980 it has been also defined [10][11] the notion of secure digital signatures.
The definition says a signature scheme is secure against a chosen-message attack by an exixtential forger if an adversary who can get valid signatures for messages of his choice m_i, is not able to produce a valid signature for any message that is different from all m_i.

## ElGamal encryption

Let G be a subgroup of prime order q of the group of prime order p and q|(p-1) (see Schnorr groups in number_theory.md).
Let g be a fixed element in G. Values p, q, g are publicly known.
A private key is randomly chosen z from Z_q, public key is e = g^z.

To encrypt m, Bob first randomly chooses r from Z_q and computes u = g^r and w = e^r * m.
Bob then sends (u, w) = (g^r, g^(r * z) * m) to Alice who deciphers it by dividing w by g^(r * z).

### Is naive ElGamal IND-CPA secure?

In chosen-plaintext attack, the adversary sends two message m0, m1 to the oracle, which chooses one, encrypts it and return to the adversary.
Now, adversary needs to be able to distiniguish which message has been encrypted with some probability 0.5 + epsilon, where epsilon is not negligible.

The encryptions of these two message would look like (g^r, g^(r * z) * m0) and (g^r, g^(r * z) * m1).

On the first look we could argue that r is uniformly random, thus g^(r * z) is uniformly random. 
How could then we distinguish g^(r * z) * m0 from g^(r * z) * m1?
The answer is it cannot be distinguished. But the problem is we have (g^r, g^(r * z) * m), not g^(r * z) * m. We need to prove IND-CPA differently.

Let's say we have an algorithm which breaks ElGamal IND-CPA. We can show that in this case we can break DDH problem (see number_theory.md).

DDH: having g, g^a, g^b, Z, check whether Z = g^(a * b).

We define encryption as in ElGamal with a as private key and (g, g^a) as public key, but instead of (g^r, g^(r * a) * m) we use (g^r, Z * m).

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/elgamal_security_proof.png

Now we encrypt m by using b as r: (g^b, Z * m). 
If Z = g^(a * b), we should be able with our algorithm to distinguish between (g^b, Z * m0) and (g^b, Z * m1).
If Z is random, our algorithm guesses the right m only with probability 0.5.

That means if there exists algorithm to break ElGamal CPA, there exists also an algorithm to break DDH. But we believe an algorithm to break DDH does not exist, so we believe algorithm to break ElGamal does not exist as well.

Note that if we have m0, m1 and get back from the encryption oracle (u, v) = (g^r, g^(r * z) * m) where m = m0 or m = m1. 
If we are able to say whether it is m0 or m1, that means we can easily compute g^(r * z) (we simply multiply v by m0_inv or m1_inv) and we solved CDH (see number_theory.md) for (g, g^r, g^z) where g^r was known only after the encryption oracle was called. But for a given (g, g^a, g^b), we can only solve DDH.

### Is naive ElGamal IND-CCA secure?

However, the naive version of ElGamal is vulnerable to chosen-ciphertext attack. 
If the adversary has a ciphertext (u, w) for some plaintext m, he can send (u, w * m1) to the decryption oracle and gets (w * m1)/u^z = m * m1 from where he can obtain m (by multiplying with m1_inv).

This can be solved by using hash function H. 
In this case the ciphertext looks (u, w, H(m)) and after decryption Alice checks whether the hash is correct (whether the encryptor knew the plaintext). 
However, this fails when the set of all possible messages is small.

As we saw above even if the adversary knows that the ciphertext has been produced by one of the two messages m_0, m_1, he still must not be able to distinguish which one was encrypted (indistinguishability-secure from chosen-ciphertext attack).

## Cramer-Shoup encryption

Cramer-Shoup encryption is indistinguishability-secure from chosen-ciphertext attack if DDH is hard (see number_theory.md) in group G and hash function H is collision-resistant.

This is called reductionist security claim - something is secure as long as some known hard problem cannot be cracked.

Previous public key cryptosystems that were secure against adaptive chosen ciphertext attack were only obtained by non-interactive zero-knowledge techniques which are much less efficient. 
On the other hand there were known practical cryptosystems but provably secure only in random oracle.

We have a group G of prime order q in the field of p elements and q|(p-1). Two random non-identity elements are publicly known.
We have a collision-resistant hash function H which takes triples of integers mod p and outputs an integer between 0 and q-1 (for example it maps 3072-bit numbers to 160-bit numbers).

We use the following notation - let x = (x1, x2), y = (y1, y2), z = (z1, z2) be pairs of integers from Z_q, let g = (g1, g2), we then say g^x = g1^x1 * g2^x2 and g^(r * x) = g1^(r * x1) * g2^(r * x2).

Alice private key consists of three randomly generated pairs x, y, z and her public key consists of the three three group elements c = g^x, d = g^y, e = g^z.

To send a message m from G, Bob chooses a random r and set:

```
u1 = g1^r
u2 = g2^r
w = e^r * m
```

Let's denote u = (u1, u2).

Bob then computes `h = H(u1, u2, w) = H(g1^r, g2^r, g1^(r * z1) * g2^(r * z2) * m)` of the concatenation of these three integers.
He computes v = c^r * d^(r * h) and sends Alice the ciphertext `(u1, u2, w, v) = (g1^r, g2^r, g1^(r * z1) * g2^(r * z2) * m, g1^(r * x1) * g2^(r * x2) * g1^(r * h * y1) * g2^(r * h * y2))`.

To decipher, Alice first computes h = H(u1, u2, w) and then uses her secret key to find u^(x + h * y) which should be equal to v. If it is not, Alice rejects the message.

If it is equal, Alice knows that Bob enciphered the message properly. This eliminates the chosen ciphertext attacks as described above for ElGamal.
Due to hash used, a valid ciphertext cannot be prepared for a message that is not known (in the attack on the naive ElGamal a valid ciphertext can be prepared for m * m1 although m * m1 is not known - the ciphertext is (g^r, g^(r * z) * m * m1, where (g^r, g^(r * z) * m is known and we just multiply the second value by m1).
Due to r used in v, the number of possible values for v is large (if only m would be hashed and a set of possible messages would be small, the system would not work).

If this test is passed, Alice decrypts the ciphertext by dividing w by u^z.
That means w/u^z = (g1^(r * z1) * g2^(r * z2) * m) / (g1^(r * z1) * g2^(r * z2)) = m.

By definition the tuple (u1, u2, w, v) is not valid if log_g1(u1) != log_g2(u2).
It can be shown that the first test will fail if the tuple is not valid (except for the very rare cases).

### Security of Cramer-Shoup

(note: the reasoning here is only a sketch of a proof) 

We will show that if we can break Cramer-Shoup, we can break DDH.

Let's say we can break Cramer-Shoup. We are given a tuple (g1, g2, u1 = g1^r1, u2 = g2^r2) and we wants to determine whether r1 = r2.

We build a simulator which encrypts as Cramer-Shoup, but uses u1 and u2 from the given tuple. 

#### INC-CPA

When we send two messages m_0 and m_1 to the simulator, it chooses bit k (which message to encrypt), encrypts m_k and returns the ciphertext to us.

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/cramer_shoup_security_proof.png

 1. If r1 = r2 the encryption is as in Cramer-Shoup and (as we are able to break it) we will be able to determine which message has been encrypted with probability which is (for some non-negligible amount) bigger than 0.5.

 2. If r1 != r2 we need to show that the returned ciphertext is independent from k.

If we prove these two lemmas, we have a statistical test to determine whether r1 = r2, meaning we broke DDH.

1 is pretty straight-forward. For 2:

Let's say g2 = g1^k.
Let's say r2 = r1 + l (we don't need to know the actual values k and l).
Message m is then masked as: 

```
g1^(r1 * z1) * g2^(r2 * z2) * m = g1^(r1 * z1 + k * r2 * z2) * m = g1^(r1 * z1 + k * r1 * z2 + k * l * z2) * m = g1^(r1 * z1) * g2^(r1 * z2) * g2^(l * z2) * m
```

While we can determine between g1^(r1 * z1) * g2^(r1 * z2) * m_0 and g1^(r1 * z1) * g2^(r1 * z2) * m_1, we cannot determine these two values multiplied by g2^(l * z2).


#### INC-CCA

Here we are allowed to send arbitrary ciphertexts (except for the challenged ciphertext) and get plaintexts.

However, we cannot send a valid ciphertext which will help us.

If we want to send an encryption of m * m1 (as in naive ElGamal), we will be able to compute hash h1 of (g1^z1 * g2^z2)^r * m * m1, but we won't be able to compute (g1^y1 * g2^y2)^(r * h1).




We need to show that we cannot construct a valid ciphertext (that will not be rejected by decryption) out of the existing ciphertext.

If we construct a ciphertext using r1, r2 (r1 != r2), it will be rejected.
To provide a valid ciphertext we would need to construct g1^(r1 * x1) * g2^(r2 * x2) * g1^(y1 * r1 * h) * g2^(y2 * r2 * h) out of knowing r1, r2 (r1 and r2 are different from the one in simulator, we choose them), g_1 ^ x1 * g_2 ^ x2, g_1 ^ y1 * g_2 ^ y2, g_1 ^ z1 * g_2 ^ z2, but we cannot construct this.

## RSA signatures 

There is a public key (n, e) and a private key (n, d), n is a product of two large primes, d is an inverse of e modulo phi(n). H is a full-domain (values range through the full interval [0, n)) hash function.

To sign a message m: `s = H(m)^d (mod n)`.
To verify a signature - the following should hold: `s^e (mod n) = H(m)`.

Now let's see that if there exists an existential forger who can successfully execute an adaptive chosen message attack in the random oracle model, the e-th root modulo n can be computed (which is believed to be hard).

For signatures, adaptive chosen message attack in random model means that a forger can send messages m_i to the Alice (forger knows only Alice's (n, e)) and when m_i is called for the first time, he obtains H(m_i) (forger cannot compute hash since it is random oracle model), while when m_i is called for the second time, he obtains signature of m_i: s_i. For some message m_i for which the signature was not obtained (only the hash), the forger can provide a valid signature.

Now we have y for which we would like to get x such that: `y = x^e (mod n)`.
We use forger as a black box - it is sending messages m_i to us. When we receive m_i for the first time we send back h_i = x_i^e (mod n) where x_i is randomly generated from Z_n. Just in one case, let's say for m_y, we send back y.
When m_i is sent for the second time, we return x_i (which a forger will find as a valid signature of m_i). 

If m_y is sent for the second time, the forger will fail. However, if it happens that m_y is the message, for which a forger provided a valid signature, we have our x. If not, we repeat the process until this happens.

## Schnorr signatures

For a longer description of Schnorr signatures see proofs1.md. However, the signer has a secret key s from Z_q, the public key is t = g^s (mod p). 

For signing, signer chooses random r from Z_q, computes x = g^r (mod p), then concatenates message m with x and computes hash (h = H(m, x)), then computes y = (r + s * h) (mod q). The signature is (y, h).

The reductionist security argument is that forging Schnorr signature is the same as breaking discrete logarithm. Why? Because the signer could any other hash function (let's say H1; to be more formal here a random oracle should be used for getting hash), compute h1 = H1(m, x) and y1 = (r + s * h1) (mod q), and finally compute (y1 - y)/(h1-h) = s * (h1 - h) / (h1 - h) = s, which is a discrete logarithm of t.

However, this security result is weak, since here is nothing about chosen-message attacks.

### Forking lemma and chosen-message existential forger for Schnorr signatures

If there exists a probabilistic (supplied with a long sequence of random bits that are used to make choices at various points) chosen-message existential forger for Schnorr signatures, then we can break discrete logarithm.

Forger is a computer program which have a non-negligible probability to produce a forged signature.

We give public key to the forger, then it can make a bounded number of hash function inquiries and signing inquiries. Meaning, for a message m_i, it gets H(m_i, x_i) for hash queries and signatures (y, H(m_i, x_i)) for signing queries.

Note that H gives a random value (random oracle) and x_i is something that is in Schnorr signing obtained as x_i = g^r_i (mod p) for some random r_i.

Now we make a guess and say that the forger will provide a valid signature at j-th inquiry. We start two instances (forking) of the forger. For the first 1, 2,..., j-1 inquiries, we answer the same for both of them. For the j-th inquiry we answer different hashes. It can be shown that there is a non-negligible probability, that both forgers produce a valid signature for the j-th inquiry. In this case we have:

```
(y = r + s * h, h)
(y' = r + s * h', h')
```

and we can extract s as (y-y')/(h-h').


[1] KOBLITZ, N. AND MENEZES, A. 2004. Another look at “provable security.” Cryptology ePrint Archive, Report 2004/152. To Appear in Journal of Cryptology. http://eprint.iacr.org/2004/152/.

[2] Oded Goldreich. On post-modern cryptography. IACR Cryptology ePrint Archive, 2006:461, 2006.

[3] S. Goldwasser, S. Micali and R. Rivest, A “paradoxical” solution to the signature problem, Proc. 25th Annual Symp. Foundations of Comp. Sci.,1984, pp. 441-448.

[4] D. Bleichenbacher, A chosen ciphertext attack against protocols based on the RSA encryption standard PKCS #1, Advances in Cryptology – Crypto ’98, LNCS 1462, Springer-Verlag, 1998, pp. 1-12.

[5] S. Goldwasser and S. Micali, Probabilistic encryption and how to play mental poker keeping secret all partial information, Proc. 14th Annual Symp. Theory of Computing, ACM, 1982, pp. 365-377.

[6] S. Goldwasser and S. Micali, Probabilistic encryption, J. Computer and System Science, 28 (1984), pp. 270-299.

[7]  C. Rackoff and D. Simon. Noninteractive zero-knowledge proof of knowledge
and chosen ciphertext attack. In Advances in Cryptology–Crypto
’91, pages 433–444, 1991.

[8] Y. Watanabe, J. Shikata and H. Imai, Equivalence between semantic security and indistinguishability against chosen ciphertext attacks, Public Key Cryptography – PKC 2003, LNCS 2567, Springer-Verlag, 2003, pp. 71-84.

[9] O. Goldreich, Foundations of Cryptography, Vol. 2, Cambridge University
Press, 2004.

[10] S. Goldwasser, S. Micali and R. Rivest, A “paradoxical” solution to the signature problem, Proc. 25th Annual Symp. Foundations of Comp. Sci., 1984, pp. 441-448. 

[11] S. Goldwasser, S. Micali and R. Rivest, A digital signature scheme secure against adaptive chosen-m
