Public-key cryptography started with Merkle puzzles (1974), Diffie-Hellman (1976), and RSA (1977).

# Diffie-Hellman

Diffie-Hellman is the first practical key exchange algorithm (it enables key exchange without an online trusted third party). The algorithm implementation is covered in Matasano challenge 33.

Fix a large prime p (600 digits or 2048 bits should be fine) and fix an integer g in {2,...,p-2}.
See Matasano challenge 35 to see why g should not be 1, p or p-1.

The communication goes like:

* Alice chooses a random a in {1,...,p-1} and then computes A = g^a (mod p)
 * Alice sends A to Bob
 * Bob chooses random b in {1,...,p-1} and computes B = g^b (mod p)
 * now the shared secret key is g^(a*b) (mod p)

Note that a and b need to be "sufficiently" large and randomly generated.
Alice can compute it as B^a (mod p), Bob as A^b (mod p).

The eavesdropper can see p, g, A = g^a(mod p), B = g^b (mod p), but he cannot compute g^(a*b) (mod p) (it would take a way too long).

## How to choose parameter g

Let us say that g is an element generating a subgroup G of Z_p* of order q. 
If g is a generator for entire Z_p*, then q = p-1, otherwise it holds: q | p-1 (Lagrange theorem).

Integer g does not need to be a generator for entire Z_p* (usually it is not), however it is crucial that q is "sufficiently" large - if subgroup G is small, the attacker could do a brute force search (the attacker needs to check all possible values which are not many if G is small from {g, g^2, g^3,...,g^(q-1)} to extract a from g^a).

The subgroup G being "sufficiently large" is actually not enough. If the order of G is not divisible by some "large" prime, then Pohlig-Hellman attack can be used (see [dlog_attacks.md](https://github.com/miha-stopar/crypto-notes/blob/master/dlog_attacks.md)).

# Public key encryption

Public key encryption is a triple of algorithms (G, E, D), where: 

 * G is a randomized algorithm that outputs a key pair (pk, sk)
 * E(pk, m) is a randomized algorithm that takes m from M and outputs c from C
 * D is a deterministic algorithm that takes c from C and outputs m from M or ⊥ (bottom symbol when authentication tag is not valid)
 * `for all (pk, sk) and all m from M: D(sk, E(pk, m)) = m`

# Semantic security

Challenger generates public/secret key and gives the public key to the adversary.
Adversary sends two messages m0, m1 (the same length) to the challenger, challenger returns encryption of one of these two messages c = E(pk, mb). 

Informally (this is not a formal definition), adversary should not be able to guess which message has been encrypted - the probability of the correct guess should be very close 0.5.

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/semantic_security_public_key.png

Note that in public key encryption there is no need to give an attacker the ability to mount the chosen-plaintext attack because an attacker has a public key and can encrypt any message (in symmetric encryption an attacker needs help from the challenger to create ciphertexts). See [security_definitions.md](https://github.com/miha-stopar/crypto-notes/blob/master/security_definitions.md).

# Chosen ciphertext security (IND-CCA - chosen ciphertext indistinguishability) (also against active attacks, not only eavesdropping)

Even if the attacker has the ability to decrypt all the ciphertexts (other than the challenged ciphertext), he cannnot decrypt the challenged ciphertext.

Challenger generates public/secret key and gives the public key to the adversary.
Adversary sends a bunch of ciphertexts c_i and the challenger returns the decrypted messages m_i = D(k, c_i). This is CCA phase 1.

Then the adversary submits two equal length messages m_0, m_1 and receives the ciphertext c = E(pk, m_b) of the first or second message.

Now the adversary can continue with the ciphertext query. This is CCA phase 2. The adversary retrieves new bunch of m_i = D(k, c_i). Obviously, c_i should not be c.

Informally, adversary should not be able to guess which message (m_0 or m_1) has been encrypted.

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/cca_security_public_key.png

# Trapdoor functions (TDF)

For public key cryptography it is crucial to have the algorithms that are easy to process in one direction, but difficult in the other (these are trapdoor functions). The bigger the difference between the difficulty of going one direction and the other, the more secure a crypto system is.

A trapdoor function X -> Y is a triple of efficient algorithms (G, F, F^(-1)).

 * G():  outputs randomized key pair (pk, sk)
 * F(pk, x): deterministic algorithm that defines a function X -> Y
 * F^(-1)(sk, y): function Y -> X that inverts F(pk, x)

For every (pk, sk) output by G and every x from X: F^(-1)(sk, F(pk, x)) = x.

Trapdoor function (G, F, F^(-1)) is secure if F(pk, x) is a "one-way" function - it can be evaluated, but cannot be inverted without sk (the probability that the adversary can guess the x is negligible).

# Public key encryption from TDFs

Three things are needed:
 * (G, F, F^(-1)): secure TDF X -> Y
 * (E_S, D_S): symmetric authenticated encryption defined over (K, M, C)
 * H: X -> K a hash function

Public key encryption system (G, E, D) is defined as:
 * key generation G as in TDF
 * E(pk, m) is defined as: take random x from X, calculate y = F(pk, x), calculate k = H(x), calculate c = E_S(k, m), outputs (y, c)
 * D(sk, (y, c)) is defined as: calculate x = F^(-1)(sk, y), calculate k = H(x), calculate m = D_s(k, c), outputs m

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/trapdoor_encrypt.png

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/trapdoor_decrypt.png

So TDF is used only to encrypt a random x, the actual message is encrypted using a symmetric encryption.

There is a theorem that says: if (G, F, F^(-1)) is a secure TDF, (E_s, D_s) provides authenticated encryption and H: X -> K is a random oracle (random function from X to K), then (G, E, D) is CCA secure. 

# RSA

RSA works due to the fact that it is easy to multiply two prime numbers, but hard to factor the product into the two primes.

## RSA trapdoor permutation

The RSA arithmetic is done in Z_N = {0,1,2,...,N-1}, where N = p * q and p, q are primes.

Remember: 

 * Z_N* is a set of all invertible elements in Z_N
 * x from Z_N is invertible if and only if gcd(x, N) = 1
 * number of elements in Z_N* is given by Euler totient function phi
 * phi(N) = (p-1) * (q-1) = N - p - q + 1

RSA trapdoor triple (G, F, F^(-1)):

### G()

G() chooses two random primes p, q and set N = p * q. The primes p and q should be selected so that factoring of n is computationally infeasible. Both need to be sufficiently large and about the same bitlength (1024 bits should be fine). There are some further requirements like p - q should not be too small.

Next, an integer e is chosen such that gcd(e, phi(N)) = 1. Lastly, an integer d is chosen such that e * d = 1 mod phi(N). Public key is a pair (N, e), private key is a pair (N, d).

Private key d in calculated using extended Euclidean algorithm (see[number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md)). For a given a and b the extended Euclidean algorithm finds x and y such that a * x + b * y = gcd(a, b). 

For RSA we want to find x and y such that phi(N) * x + e * y = gcd(phi(N), e). An integer e was chosen such that gcd(phi(N), e) = 1, thus there exists x and y such that phi(N) * x + e * y = 1. Apply mod phi(N) and you get e * y = 1. Number y is what we are looking for: private key d.

Output pk = (N, e), sk = (N, d)

### F(pk, x)

`F(pk, x) = x**e % N`

Message x should be represented as an integer from [0, n-1].

### F^(-1)(sk, y)

`F^(-1)(sk, y) = y^d % N`

It holds: `F^(-1)(sk, y) = y^d = x^(e*d) = x^(k*phi(N)+1) = (x^phi(N))^k * x = x`

Euler theorem was used in the last step: for every x from Z_N* it holds that `x^phi(N) = 1`. Obviously, the simple proof above is invalid when x not from Z_N*. However, the equation still holds, see below:

If x and N are not coprimes (but it is very difficult to find an integer that is not relatively prime to N), then `x = 0 mod p` or `x = 0 mod q` (due to N = p * q). Due to symmetry we can assume x = 0 mod p.

We can write as above:

`x^(e*d) = x * x ^ (k * phi(N)) = x * x ^ (k * (q-1) * (p-1)) = x * (x^(q-1))^(k*(p-1) = x mod q`

The last step follows from Fermat little theorem: `x^q = x mod q` where q is prime. If you divide this equation by x you get the desired property: `x ^ (q-1) = 1 mod q`.

Now we know:

 * `x ^ (e*d) = x mod q`
 * `x ^ (e*d) = 0 mod p` (because x = 0 mod p)

This means `x ^ (e*d) - x` is a multiple of both p and q. As N = p * q (and p and q are primes) we see that `x ^ (e*d) - x` is a multiple of N as well. Thus, `x ^ (e*d) = x mod N`.

NOTE: you should never ever directly encrypt with RSA! You should use it only as a trapdoor function.

Padding, like PKCS1 or OAEP, is applied to randomize the message (and to prevent message to be too small - m**e must not be less than N, otherwise m can be reverted efficiently).

## RSA problems

### Factoring attack

Never should be used the same RSA modulus (N) in two different public keys. Even more, all RSA modulus must be pairwise coprime (relatively prime), otherwise the factorization of modulus get much easier. See below.

Let us say we have N1 = p*q1 and N2 = p*q2. Now, we can calculate gcd(N1, N2) using Extended Euclid Algorithm (see [number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md)) which is fast (as opposed to the factorization).

We will get gcd(N1, N2) = p. Now we can trivially calculate q1 and q2:

 * `q1 = N1/p`
 * `q2 = N2/p`

Thus, we have `phi(N1) = (p-1) * (q1-1)` and `phi(N2) = (p-1) * (q2-1)`. Using phi(N1) and phi(N2) we can calculate private keys (N, d1), (N, d2) for both public keys (N, e1), (N, e2) by calculating inverse of e1 and e2 modulo phi(N1) and phi(N2) respectively.

### Chosen-ciphertext attack on plain RSA

Plain RSA (without PKCS1 or OAEP padding) is susceptible to a chosen-ciphertext attack. It is trivial:

 * attacker wishes to decrypt c which was encrypted c = m^e (mod N) 
 * it holds `m = c^d (mod N)`
 * attacker chooses a random integer s
 * attacker calculates: `c1 = s^e * c`
 * attacker asks for a decryption of a message c1
 * attacker gets a decrypted c1: `m1 = c1^d (mod N)`
 * `m1 = c1^d (mod N) = (s^e * c)^d (mod N) = s^(e*d) * c^d (mod N) = s * c^d (mod N) = s * m (mod N)`
 * attacker calculates: `m = m1 * s_inv (mod N)` where s_inv is inverse of s modulo N

# ElGamal

Public key encryption can be based on trapdoor functions (such as RSA) or Diffie-Hellman protocol. ElGamal is built on Diffie-Hellman.

Diffie-Hellman can be converted into public key encryption as follows:

 * Alice chooses a random number a from {1,...,n}
 * h=g^a is treated as a public key, a is treated as a secret key
 * Bob computes g^(ab) (he chosed number b from {1,...,n})
 * Bob derives key k from g^(ab) and encrypts message m with k

Formally:

 * G: finite cyclic group of order n
 * (E, D): symmetric auth. encryption defined over (K, M, C)
 * H: GxG -> K is a hash function

Note: h = g^a

E(pk=(g, h), m):

 * choose b from Z_n, calculate u = g^b, calculate v = h^b
 * calculate k = H(u, v) and c = E(k, m)
 * output (u, c)

Note: g^b needs to be publicly known - it is needed for decryption

D(sk=a, (u, c)):

 * calculate v = u^a
 * calculate k = H(u, v) and m = D(k, c)
 * output m


# DSA

## Key generation

A prime number q is chosen so that: 2^159 < q < 2^160. This means q has bit length 160. A prime number p is chosen so that: 

 * 2^(511+64j) < p < 2^(512+64j) for some j from {0,1,...,8}
 * q | p - 1

The bit length of p is a multiple of 64 and is between 512 and 1024.

Now we select a generator g of a subgroup of order q in Z_p\*. For each divisor of p-1 there exists a subgroup of this order (see Cauchy in [number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md)), for each subgroup there is at least one generator. The generator can be found by the procedure described in [number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md): choose x from Z_p\* such that x\*\*((p-1)/q) is not 1 in Z_p.

Now select a random integer a from {0, 1, 2,...,q-1} and compute A = g\*\*a mod p.

The public key of a signer is (p, q, g, A), the private key is a. The number q is about 2\*\*160 which means that computing the secret key a from A requires solving the discrete logarithm problem in the subgroup of order about 2\*\*160.

## Signature generation

A binary message x of an arbitrary length can be signed. First a hash function is applied, like SHA-1:

`SHA-1: {0,1}* -> {0,1}**160`

The signer then chooses a random number k from {1,2,...,q-1} and computes:

`r = (g**k mod p) mod q`

Then the signer sets:

`s = (k**(-1) * (SHA-1(x) + a * r)) mod q`

The expression k\*\*(-1) means the inverse of k modulo q. The inverse modulo q can be calculated using Extended Euclidean algorithm, see [number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md).

The signature of x is (r, s).

## Verification

The verifier first verifies that 1 <= r <= q-1 and 1 <= s <= q-1. If not, the verifier rejects the signature. Let us denote s**(-1) as s_inv.

It holds:
```
k = s_inv * (SHA-1(x) + a*r)
g**(s_inv * (SHA-1(x) + r*a)) = g**(s_inv * SHA-1(x)) * g**(s_inv * r * a) = g**((s_inv * SHA-1(x)) mod q) * A**((s_inv * r) mod q)
```
Thus, the verifier needs to verify that:

`r = ((g**((s_inv * SHA-1(x)) mod q) * A**((r * s_inv) mod q)) mod p) mod q`

## Matasano challenges related to DSA

### Matasano challenge 43

The random number k used when signing is sometimes refered as nonce, but is actually not a nonce in a traditional sense (the value that never repeats), because in DSA this number must also be unpredictable and kept secret.

If the number k is known, it is trivial to extract the private key:

```
s = (k_inv * (SHA-1(x) + a * r)) mod q
s * k = (SHA-1(x) + a * r) mod q
a = ((s * k - SHA-1(x)) * r_inv) mod q
```

Also, if the k is chosen from a small pool of numbers, all these numbers can be tried - you just take a candidate for k, calculate private key a, sign a message, and compare the signature that is known and signature that you just calculated. If the same, you know k and private key a.

### Matasano challenge 44

Even if k is unknown, it can be easily computed if two messages are signed with the same k and the same private key:

```
s_1 = k_inv * (SHA-1(x_1) + a * r)
s_2 = k_inv * (SHA-1(x_2) + a * r)
```

The same private key a is used for signing messges x_1 and x_2. The value r is the same as it depends only on k. If we multiply both equations with k, we get (modulos omitted):

```
s_1 * k = (SHA-1(x_1) + a * r)
s_2 * k = (SHA-1(x_2) + a * r)
```

Now we use subtraction:

`s_1 * k - s_2 * k = (SHA-1(x_1) + a * r) - (SHA-1(x_2) + a * r)`

And from there:

```
k * (s_1 - s_2) = SHA-1(x_1) - SHA-2(x_2)
k = (s_1 - s_2) / (SHA-1(x_1) - SHA-2(x_2))

```

Note that once again when dividing you have to find an inverse modulo q of the denominator and then multiply this value with the numerator.

### Matasano challenge 45

Parameters p, q, and g must be properly chosen. Being g for example 0 modulo p is wrong, as this would mean that public key A is 0. Similarly, g being 1 modulo p would mean A is 1. This is because:

```
A = (g**a) % p = (g * g * g * ... * g) % p = (g % p) * (g % p) * (g % p) * ... * (g % p) = 1 * 1 * 1 * ... * 1 = 1
```

In Matasano challenge 45 the scenario where the verifier can be forced to set parameter g to the value p+1 is discussed. Anybody can forge any message to have a valid signaturethen. Let A be the public key of the signer and h the hash of some message.

```
r = ((A ** h) % p) % q
s = (r * h_inv) % q
```

The pair r, s is a valid signature for a message with hash value h. Thus, the valid signature has been calculated ONLY using the public key. This is because in the equation that needs to be checked, the first factor is 1 (as ((p+1)**anything) % p = 1):

`r = ((g**((s_inv * SHA-1(x)) mod q) * A**((r * s_inv) mod q)) mod p) mod q`


Thus, what needs to be checked is: 

`r = ((A**((r * s_inv) mod q)) mod p) mod q`

Our r is:

`r = ((A ** h) % p) % q`

Thus we need to check whether `A**h = A**((r * s_inv) mod q)`. Or the same - whether `h = (r * s_inv) mod q`.

And this equation holds because of how we defined s.






