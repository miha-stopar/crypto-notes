# Key Derivation Function (KDF)

How to derive many keys from one? Derivation of many keys from one source key is needed in practice all the time and is done using Key Derivation Functions (KDFs).

Usually a source key (SK) is sampled from hardware random number generator or key exchange protocol.

When source key is uniform a PRF can be used to construct KDF:
KDF(SK, CTX, L) := F(SK, (CTX || 0)) || F(SK, (CTX || 1)) || ... || F(SK, (CTX || L))

where F is a secure PRF, CTX is a string that uniquely identifies the application (to prevent different applications to have the same keys), L is length.

But what if the source key is not uniformly distributed in K space? In this case PRF output may not look random. And unfortunately the keys are often not uniformly random (by key exchange protocol the key is uniform in just some subset of K; hardware RNG may produce biased output).

The solution is extract-then-expand paradigm: extract pseudo-random key k from source key (distribution indistinguishable from uniform) and then expand k using PRF as above.

An extractor usually takes as input a salt - fixed non-secret string chosen at random.

## HKDF

HKDF is a KDF based on HMAC (HMAC is used both as PRF and as an extractor). It implements the extract-then-expand paradigm.

 1) k = HMAC(salt, SK)
 2) expand using HMAC as a PRF with key k

In the first step salt is used as HMAC key and SK as HMAC data. In the second HMAC is used as a PRF described above to generate a key k (with as many bits as needed).

Note: once you obtain source key (either from hardware or from key exchange protocol), you never use it directly as a session key - you must run SK through the KDF and thus get the keys to be used.

# Password-Based KDF (PBKDF)

Passwords have generally low entropy - not enough to generate session keys.

Do not use HKDF to generate session keys out of passwords as this is easily compromised by dictionary attacks. The standard approach is PKCS#5 where a salt and a slow hash function are used (key is guessable, but guessing is very slow due to slow hash function):
H_(c)(pwd || salt): iterate hash function c times

This is why you should never use general purpose hash functions (like SHA family) for deriving a key out of password - these hash functions are fast and the attacker can use dictionary attacks or brute force attacks. You should use something like bcrypt or scrypt instead.

# Merkle-Damgard paradigm

From a given collision resistant function for short messages, how to construct collision resistant (collision resistance is described in [security_definitions.md](https://github.com/miha-stopar/crypto-notes/blob/master/security_definitions.md)) function for long messages? All standard collision resistant hash functions follow Merkle-Damgard paradigm to achieve this.

In Merkle-Damgard iterated construction we assume we have h: T x X -> T where h is collision resistant hash function for small inputs. 
Function h is called one-way compression function because it transforms two fixed-length inputs into a fixed-length output and it is difficult to compute inputs which compress to a given output.

First, the message m is split into blocks m[0],...,m[k]. The last block is padded: m[k] || PB. 
Padding block is defined as 100...0 || msg len (some fixed number of bits is reserved for message length, see for example Expandable messages attack in attacks.md to see why appending message length is important). If there is no space for PB another block is added.
There is a fixed (always the same) IV.

```
H_0 = IV
H_1 = h(H_0, m[0])
H_2 = h(H_1, m[1])
...
H_k = h(H_(k-1), m[k-1])
H_(k+1) = h(H_k, m[k] || PB)
```

The variables H_0, H_1,... are called chaining variables.

It can be proved that if h is collision resistant then so is H (it can be quickly shown that if a collision exists for H, there is a collision for h as well - thus, no collision exists for H which is not a consequence of collision of h).

## Davies-Meyer

Compression functions are often built from block ciphers (like SHACAL), an example is Davies-Meyer compression function.

If (E, D) is a block cipher over (K, X) where X = {0,1}^n. Davies-Meyer compression function derived from E is defined as:

`h(x, y) = E(y, x) xor x`

Using Davies-Meyer in Merkle-Damgard paradigm in each step the inputs are a chaining variable H_i and message block m[i].

`h(H_i, m[i]) = E(m[i], H_i) xor H_i`

The message block is used as a block cipher key. Usually, block ciphers are optimized to encrypt long messages with a fixed key - changing the key on every block might be very slow. Thus, not all block ciphers can be used with Davies-Meyer (use for example SHACAL instead of AES).

AES should not be used in Davies-Meyer also from other reasons. Using AES would mean the output is 128-bit long, which is not enough for collision resistance - a collision could be found with 2^(64) evaluations, see birthday paradox in [security_definitions.md](https://github.com/miha-stopar/crypto-notes/blob/master/security_definitions.md).

Also, it is more efficient to use block ciphers which use longer keys (at least 512, unlike AES) because this way more message bits are processed in every round.

For example SHA-256 is built on Merkle-Damgard paradigm, using Davies-Meyer compression function with SHACAL-2 block cipher.





