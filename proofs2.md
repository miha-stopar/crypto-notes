# Shamir's secret sharing scheme

A secret is divided into n parts - one, unique part for each participant. For a threshold k < n it holds: from any k parts a secret can be reconstructed (k-1 or fewer parts leave secret completely undetermined).

Let's say a_0 is a secret. We then randomly choose k-1 coefficients: a_1, a_2, ..., a_(k-1).
We now have a polynomial of degree k-1:

```
f(x) = a_0 + a_1 * x + a_2 * x^2 + ... + a_(k-1) * x^(k-1)
```

The first coefficient a_0 is our secret. Now we compute n points, let's say f(x_1), f(x_2), ..., f(x_n).
We distribute points to the participants (one point for one participant). 
Now if any of the k points are put together, f(x) can be computed because k points determine the polynomial of a degree k-1 (you have k equations with k variables). Lagrange interpolation can be used to compute f. Once the f is known, secret is computed as f(0).

However, if natural numbers are used for a_i, then already k-1 points reveal some information - you have k-1 equations and k variables, which means you can choose for one variable whatever you want and you will still be able to find a solution for k-1 equations. However, a solution needs to be in natural numbers which limits all possible choices.
To avoid this we use coefficients from a finite field Z_p where p is prime bigger than all coefficients:

```
f(x) = a_0 + a_1 * x + a_2 * x^2 + ... + a_(k-1) * x^(k-1) (mod p)
```

## Verification of shares

Each participant who received (x_i, f(x_i)) can verify whether the share is correct.
Points x_1, ..., x_n are publicly known, as well as are A_0 = g^a_0 (mod p), A_1 = g^a_1 (mod p), ..., A_(k-1) = g^(k-1) (mod p).

When i-th participant receives (x_i, f(x_i)), she checks whether:

```
g^(f(x_i)) = g^(a_0 * a_1 * x_i + ... + a_(k-1) * x_i^(k-1)) = A_0 * A_1^x_i * ... * A_(k-1) * x_i^(k-1))
```

Note that this can be checked only by the participant who knows f(x_i), third-parties cannot check whether the sharing was correct.

# Publicly verifiable encryption (not CCA secure)

Stadler [1] proposed how to build publicly verifiable Shamir's secret sharing.

How can third-party verify the sharing?
When dealer prepares shares (x_i, s_i := f(x_i)), it decrypts s_i using public key of participant P_i and sends enc_P_i(s_i) to P_i. Let's say it uses ElGamal:

```
y = h^z // z is P_i's secret key
enc_P_i(s_i) = (h^alpha, s_i_inv * y^alpha) // here inverse of s_i is used intead of s_i as normally in ElGamal, but it doesn't really matter
```

Let's make all S_i = g^s_i (mod p) publicly known.

Now dealer needs to prove to a third-party that enc_P_i(s_i) is encryption of log_g(S_i).
Somehow extended Schnorr's protocol can be used for this. In the first step, dealer does not send only g^w (mod p), but also (g^y)^w.
In the last step third-party verifier does not check only that dealer knows alpha (log_g(h^alpha)), but also checks that `(g^y)^r = (g^y)^(w + c * alpha) = (g^y)^w * ((g^y)^alpha)^c`.

Note that (g^y)^alpha is given by S_i to the power of the second element in enc_P_i(s_i): `S_i^(s_i_inv * y^alpha) = g^(s_i * s_i_inv * y^alpha) = (g^y)^alpha`.

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/stadler_verifiable_encryption.png

Basically, this is Schnorr protocol for proving the knowledge of log to base g of h^alpha and knowledge of log to base g^y of (g^y)^alpha and that in both cases the value (alpha) is the same.

The protocol can be turned into non-interactive one using Fiat-Shamir heuristics.

Note: naive ElGamal is used which is not CCA secure.

# Applications of publicly verifiable encryption

## Key escrow

Let's say we have end-to-end encryption where for each file a session key for symmetric encrytpion is randomly generated. File is encrypted using session key, session key is encrypted using party's A public key A_pub. Party A sends encrypted file and encrypted session key to the server in an untrusted domain. When party A needs a file, it downloads it from the server, decrypts session key using A_priv and then decrypts a file.

When party A wants to share a file to party B, it encrypts session key using B_pub and stores it on the server from where party B can get it and encrypts it using B_priv.

Now, what if party A (or any other party) loses its private key (A_priv)?

Each session key s_f might be encrypted using T's public key (T_pub(s_f)) and sent to T. This means that T can decrypt all files and cannot be place in an untrusted domain (it plays the role of the trusted third party now).

There were some techniques proposed how to prevent T from mass wiretapping, for example Encapsulated Key Escrow [2], however individual wiretapping is still possible.

Key Escrow where public verifiable encryption is used is for example Auto-recoverable auto-certifiable cryptosystems [3] where escrow authority have a proof that it can recover private keys.

# Paillier's encryption

We will need (modified) Paillier encryption later to build CCA2 secure verifiable encryption [4].

Paillier [6] proposed in 1999 a new public key scheme which is based on Composite Residuosity Class Problem which is related to Decision Composite Residuosity Assumption (see number_theory.md). 

Let n = p * q be the product of two large primes. In this case phi(n) = (p-1) * (q-1). By definition Carmichael's function lambda is lambda(n) = lcm(p-1, q-1).

It holds |Z_n^2*| = phi(n) * n and for any w from Z_n^2*:

```
w^lamda(n) = 1 mod n
w^(n * lamda(n)) = 1 mod n^2
```

It holds:

```
gcd(lambda(n), n) = 1
```

Consider group Z_n^2* and subgroup P of all n-th powers of elements. DCRA says it is hard to distinguish random elements of Z_n^2* from random elements of P.

Now, let's define function epsilon:

```
epsilon: Z_n x Z_n* -> Z_n^2*
epsilon: (x, y) -> g^x * y^n mod n^2
```

where g is an element of order n or some multiplication of n. Let the set of such g be B.

Firstly, how can we get an element of order n? Note that the element g = 1 + n has order n (using binomial theorem it can be quickly shown that (1 + n)^x (mod n^2) = (1 + x * n) (mod n^2) for 0 <= x < n.

Both, Z_n x Z_n\* and Z_n^2* have the same number of elements: n * phi(n). Function epsilon is actually bijective (it can be shown that it is injective).

For w from Z_n^2\*, we call n-th residousity class of w with respect to g the unique integer x from Z_n for which there exists y from Z_n\* such that:

```
epsilon(x, y) = w
```

We denote it as [[w]]_g.

Some properties of [[w]]_g:

 * [[w]]_g = 0 if and only if a is n-th residue modulo n^2
 * [[w]]_g is a homomorphism from (Z_n^2*, *) to (Z_n, +): [[w1 * w2]]_g = [[w1]]_g + [[w2]]_q mod n

The n-th Residuosity Class Problem of base g, denoted Class[n, g] is defined as the problem of computing the class function in base g: for a given w from Z_n^2* compute [[w]]_g from w.

## Class[n, g] is random-self-reducible over w from Z_n^2*

We can compute [[w']]_g. Now we want to compute [[w]]_g.

The paper [6] says it is easy to get alpha from Z_n and beta from Z_n* such that: `w' = w * g^alpha * beta^n mod n^2`.

If this is the case then:

```
w' = g^x * y^n mod n^2
w * g^alpha * beta^n = g^x * y^n mod n^2
w = g^(x-alpha) * (y * beta_inv)^n mod n^2

```

However, can we really find alpha, beta such that w' * w_inv = g^alpha * beta^n mod n^2? Isn't this [[w' * w_inv]]_g?

## Class[n, g] is random-self-reducible

Note that Class[n, g] is random-self-reducible over g from B (where B is a set of all elements of order n * k) - if we can get [[w]]_g, we can get also [[w]]_g1, because we can get(x, y) and (x1, y1) such that:


```
w = g^x * y^n mod n^2
g1 = g^x1 * y1^n mod n^2
```

Now from the second equality:

```
g^x1 = g1 * y1^(-n) mod n^2
g^x = g1^(x * x1_inv) * y1^(-n * x * x1_inv) mod n^2
```

And:

```
w = g1^(x * x1_inv) * (y1^(x * x1_inv) * y)^n mod n^2
```

where x1_inv is inverse of x1 modulo n. Note that `g^(x1 * x1_inv) = g^(k * n + 1) = g^(k * n) * g = g mod n^2` because g is of order n.

This means `[[w]]_g1 = x * x1_inv mod n`.

Thus, Class[n, g] is independant from g and we can denote it as Class[n]: given w from Z_n^2* and g from B, compute [[w]]_g.

## If we can factor n, we can get [[w]]_g (meaning Fact[n] => Class[n])

Let's see if we can compute [[w]]_(1+n) and [[g]]_(1+n) (1+n is from B)? 

There exists x1, y1, x2, y2, x3, y3 such that:

```
w = (1+n)^x1 * y1^n mod n^2
g = (1+n)^x2 * y2^n mod n^2
w = g^x3 * y3^n mod n^2
```

Let's check w^lambda(n) and g^lambda(n) (use binomial theorem and Carmichael's results mentioned above):

```
w^lambda(n) = (1+n)^(x1 * lambda(n)) * y1^(n * lambda(n)) = (1 + x1 * lambda(n) * n) mod n^2
g^lambda(n) = (1+n)^(x2 * lambda(n)) * y2^(n * lambda(n)) = (1 + x2 * lambda(n) * n) mod n^2
```

That means we can easily compute x1 and x2 if we know lambda(n):

```
x1 = (w^lambda(n) - 1)/(lambda(n) * n) mod n^2
x2 = (g^lambda(n) - 1)/(lambda(n) * n) mod n^2
```

However, we are trying to find x3. As it turns out, `x3 = x1 * x2_inv`. Why? Let's use g expressed with (1+n) in w expressed with g.

```
w = g^x3 * y3^n = ((1+n)^x2 * y2^n)^x3 * y3^n mod n^2
w = (1+n)^(x2 * x3) * (y2^x3 * y3)^n mod n^2
```

If we compare the last equation with the one where w is expressed with (1+n), we see that x1 = x2 * x3.

So we can get [[w]]_q as:

```
[[w]]_g = (w^lambda(n) - 1)/(g^lambda(n) - 1)
```

## How to get a generator

It can be shown that when gcd((g^lambda(n) - 1)/n, n) = 1, g is of order alpha * n.

```
g = (1 + n)^x * y^n mod n^2
g^lambda(n) = (1 + n)^(x * lambda(n)) * y^(n * lambda(n)) mod n^2
g^lambda(n) = (1 + n * x * lambda(n)) mod n^2
(g^lambda(n) - 1)/(n * lambda(n)) = x mod n^2

gcd(x * lambda(n), n) = gcd(x, n) // because gcd(lambda(n), n) = 1
```
Let's have a look at:

```
g^lambda(n) = (1 + n * x * lambda(n)) mod n^2
g^(lambda(n) * b) = (1 + n * x * b * lambda(n)) = 1 mod n^2
```

This means `n^2 | n * x * b * lambda(n)`. This means n | x * b * lambda(n). Because gcd(n, lambda(n)) = 1, this means n | x * b. If gcd(x, n) = 1, n | b. 

So, if gcd(x, n) = 1, then we have a generator. As we saw, x = (g^lambda(n) - 1)/(n * lambda(n)) mod n^2.

# Paillier encryption made CCA2 secure

Cramer and Shoup proposed in [5] how to make Paillier encryption CCA2 secure.

Due to fundamental theorem of finite abelian groups (see [number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md)) Z_n^2* can be decomposed as an internal direct product of cyclic groups:

```
Z_n^2* = G_n x G_n' x G_2 x T
```

where G_n is of order n (= p * q), G_n' is of order n' (= p' * q' = (p-1) * (q-1) / 4), G_2 and T are of order 2. T is generated by n^2 - 1.

This decomposition is unique except for the choice of G2 (two possible choices).

Element h = 1+n has order n and generates G_n.

Let's check `X = G_n x G_n' x T`. This is a cyclic group (the product of two cyclic groups is cyclic group when their orders are relatively prime). Let L be the subgroup of n-th powers of X.

We can see that `L = G_n' x T`. 

So L are all n-th powers in X. It can be shown that DCRA (see number_theory.md) implies it is hard to determine whether an element is from L. Cramer-Shoup modification of Paillier encryption is based on this, however below is presented Camenisch-Shoup modification of Paillier encryption, because it supports verifiable encryption and verifiable decryption.

# Publicly verifiable encryption (CCA secure)

Camenisch and Shoup [4] proposed in 2003 the first verifiable encryption system that provides CCA security. It is based on a number of techniques which will be presented first:

 * Fujisaki-Okamoto's method for proving relations on committed values (ref todo)
 * related interval proofs (ref todo)
 * Paillier encryption [6]
 * Cramer and Shoup's universal hash proof encryption technique [5]
 * Cramer, Damgard, Schoenmakers' proofs of partial knowledge (ref todo)
 * Boudot's exact interval proofs (ref todo)
 * new protocols for proving the inequality of discrete logarithms.

Questions (todo): 
why not regular hash
why 2*t in m1 := (e/u^x1)^(2*t) when decrypting

Camenisch and Shoup [4] does not fully fit into [5] and an independent security analysis was thus provided in [4].

Again, `Z_n^2* = G_n x G_n' x G_2 x T`.

The n-th power elements are `P = G_n' x G_2 x T`.

Value lambda is security parameter, l = l(lambda) is an auxiliary parameter.

H is a collision resistant hash function which maps a triple (u, e, L) to a number in the set [2^l]. The paper says hash function needs to be keyed, but seems that this is not used anywhere in the security analysis?

Let abs be a function:

```
abs: Z_n^2* -> Z_n^2*
abs: a -> (n^2 - a) if a > n^2/2
abs: a -> a if a < n^2/2
```

It holds `v^2 = (abs(v))^2' for all v from Z_n^2*.

## Key generation

Choose p, q, p', q' be distinct large primes and p = 2 * p' + 1, q = 2 * q' + 1. Primes p' and q' are of l-bit length.

Let n = p * q and n' = p' * q'.

Choose x1, x2, x3 randomly from {0, 1, 2,..., n^2/4 - 1} when n^2/4 is an integer (round down). Choose g' randomly from Z_n^2* and compute g = (g')^(2 * n), y1 = g^x1, y2 = g^x2, y3 = g^x3.

The public key is `(hashFunc, n, g, y1, y2, y3)`. The secret key is `(hashFunc, n, x1, x2, x3)`.

## Encryption

The input for encryption are plaintext m and label L. Label is specifying the context in which the encryption or decryption operation is to take place.

```
u = g^r mod n^2 // r is chosen randomly from [n/4] (smaller than n/4 - 1)
e = y1^r * h^m mod n^2
v = abs((y2 * y3^hashFunc(u, e, L))^r)
```

## Decryption

Check whether `abs(v) = v`. If not, reject.

Check whether `u^(2 * (x2 + hashFunc(u, e, L) * x3)) = v^2 mod n^2`. If not, reject.

Compute (x1_inv is inverse of x1 modulo n^2):

```
m1 = e * u^x1_inv % n^2
```

Check whether m1 is of form h^m (h is 1+n, h^m is 1+m*n). If not, reject. If yes, compute m from it.

## Camenisch-Shoup Paillier encryption is verifiable encryption

Camenisch and Shoup [4] proposed an encryption scheme which is CCA2 secure and can be combined with a proof showing that the message m (for which the decryptor knows the ciphertext) has some particular property. The property here is being discrete logarithm of some value that encryptor reveals to the decryptor.

We choose a cyclic group of prime order ro generated by gamma (both publicly known). 
We choose two generators g1, h1 of Z_n* subgroup of order n1. Note that here a different n and n1 might be used from the one in encryption, however we assume the same (the paper says it can be the same). 

We choose two additional security parameters k and k1 such that 2^(-k) and 2^(-k1) are negligible, 2^k < min{p1, q1, ro}, and ro < n * 2^(-k-k1-3).

Encryptor computes delta = gamma^m and wants to prove that (u, e, v) is encryption of log_gamma(delta).

The prover (encryptor) chooses a random s smaller than n/4, computes l1 = g1^m * h1^s and sends l1 to the verifier. The prover then chooses r1 from [-n * 2^(k+k1-2), n * 2^(k+k1-2)], s1 from [-n * 2^(k+k1-2), n * 2^(k+k1-2)], and m1 from [-ro * 2^(k+k1), ro * 2^(k+k1)]. The prover computes:

```
u1 = g^(2*r1)
e1 = y1^(2*r1) * h^(2*m1)
v1 = (y2 * y3^hash(u, e, L))^(2*r1)
delta1 = gamma^m1
l1 = g1^m1 * h1^s1
```

and sends these values to the verifier.

The verifier chooses a random challenge c < 2^k and sends c to the prover.
The prover replies with r_tilde = r1 - c * r, s_tilde = s1 - c * s, and m_tilde = m1 - c * m (computed in Z, can be negative).

The verifier checks whether:

```
u1 = u^(2*c) * g^(2*r_tilde)
e1 = e^(2*c) * y1^(2*r_tilde) * h^(2*m_tilde)
v1 = v^(2*c) * (y2 * y3^hash(u, e, L))^(2*r_tilde)
delta1 = delta^c * gamma^m_tilde
l1 = l^c * g1^m_tilde * h1^s_tilde
-n/4 < m_tilde < n/4
```

Note that this way the prover proves that he knows such r, m, s that: 

```
u^2 = g^(2*r)
e^2 = y1^(2*r) * h^(2*m)
v^2 = (y2 * y3^hash(u, e, L))^(2*r)
delta = gamma^m
l = g1^m * h1^s
-n/2 < m < n/2
```

Using Camenisch-Stadler notation [7], this can be written:

```
PK{(r, m, s): u^2 = g^(2*r) ∧ e^2 = y1^(2*r) * h^(2*m) ∧ v^2 = (y2 * y3^hash(u, e, L))^(2*r) ∧ delta = gamma^m ∧ l = g1^m * h1^s ∧ -n/2 < m < n/2}

```

Furthermore (and crucial for verifiable encryption), it can be shown that only one m satisfies these equations - (u, e, v) is really encryption of log_gamma(delta).


[1] M. Stadler, Publicly verifiable secret sharing, Advances in Cryptology — EUROCRYPT ’96, LNCS, vol. 1070, Springer Verlag, 1996, pp. 191–199.

[2] M. Bellare and S. Goldwasser, Encapsulated key escrow, Preprint, 1996

[3] A. Young and M. Young, Auto-recoverable auto-certifiable cryptosystems., Advances in Cryptology — EUROCRYPT ’98, LNCS, vol. 1403, Springer Verlag, 1998, pp. 17–31.

[4] J. Camenisch and V. Shoup, Practical verifiable encryption and decryption
of discrete logarithms, http://eprint.iacr.org/2002/161, 2002.

[5] R. Cramer and V. Shoup, Universal hash proofs and a paradigm for adaptive chosen ciphertext secure public-key encryption, http://eprint.iacr.org/ 2001/108, 2001.

[6] P. Paillier, Public-key cryptosystems based on composite residuosity classes, Advances in Cryptology — EUROCRYPT ’99 (J. Stern, ed.), LNCS, vol. 1592, Springer Verlag, 1999, pp. 223–239.

[7] J. Camenisch and M. Stadler, Efficient group signature schemes for large groups, Advances in Cryptology — CRYPTO ’97 (B. Kaliski, ed.), LNCS, vol. 1296, Springer Verlag, 1997, pp. 410–424.



