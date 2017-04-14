# Crypto notes

I started taking these notes when doing Dan Boneh crypto course on Coursera and when solving Matasano challenges. Both are highly recommended. While trying to learn some more crypto I am gradually adding further content.

A helpful resources I sometimes use:
 
 * [Handbook of applied cryptography](http://cacr.uwaterloo.ca/hac/) by Alfred J. Menezes, Paul C. van Oorschot and Scott A. Vanstone
 * [Guide to Elliptic Curve Cryptography](http://cacr.uwaterloo.ca/ecc/) by Darrel Hankerson, Alfred Menezes, and Scott Vanstone
 * [A Graduate Course in Applied Cryptography](https://crypto.stanford.edu/~dabo/cryptobook/draft_0_2.pdf) by Dan Boneh and Victor Shoup
 * [A computational Introduction to Number Theory and Algebra](http://www.shoup.net/ntb/ntb-v2.pdf) by Victor Shoup

## Content

Some basic number theory and its applications in cryptography:
[number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md).

Terms like semantic security, chosen-plaintext attack, chosen-ciphertext attack, MACs, collision resistance, random oracle model:
[security_definitions.md](https://github.com/miha-stopar/crypto-notes/blob/master/security_definitions.md).

Public key security definitions and most known schemes (ElGamal, RSA, DSA):
[public_key_security.md](https://github.com/miha-stopar/crypto-notes/blob/master/publick_key_security.md).

A brief history of provable security:
[provable_security.md](https://github.com/miha-stopar/crypto-notes/blob/master/provable_security.md).

Commitment schemes, sigma protocols, zero knowledge proofs:
[proofs1.md](https://github.com/miha-stopar/crypto-notes/blob/master/proofs1.md).

Shamir's secret sharing scheme, Paillier's encryption, publicly verifiable encryption:
[proofs2.md](https://github.com/miha-stopar/crypto-notes/blob/master/proofs2.md).

Brief descriptions of some of the attacks included in Matasano challenges (length-extension attacks, Bleichenbacher, hash collisions, Wang, RC4 biases ... ):
[attacks.md](https://github.com/miha-stopar/crypto-notes/blob/master/attacks.md).

Briefly about Elliptic Curves Cryptography (ECC) and how to speed up ECC operations:
[ecc.md](https://github.com/miha-stopar/crypto-notes/blob/master/ecc.md).

Some attacks on discrete logarithm (index calculus, Pohlig-Hellman, Pollard's rho):
[dlog_attacks.md](https://github.com/miha-stopar/crypto-notes/blob/master/dlog_attacks.md).

Briefly about pairing-based crypto:
[pairing_based_crypto.md](https://github.com/miha-stopar/crypto-notes/blob/master/pairing_based_crypto.md).

Ciphers (only AES for now):
[ciphers.md](https://github.com/miha-stopar/crypto-notes/blob/master/ciphers.md).

Some crypto constructs (key derivation functions, Merkle-Damgard paradigm, Davies-Meyer compression function):
[crypto_constructs.md](https://github.com/miha-stopar/crypto-notes/blob/master/crypto_constructs.md).

Brief description of some crypto libraries (OpenSSL, scrypt, NaCl, SJCL):
[crypto_libraries.md](https://github.com/miha-stopar/crypto-notes/blob/master/crypto_libraries.md).

How some applications use crypto:
[applications.md](https://github.com/miha-stopar/crypto-notes/blob/master/applications.md).

Some basics about trusted computing:
[trusted_computing.md](https://github.com/miha-stopar/crypto-notes/blob/master/trusted_computing.md).  

Some basic info on how to work with bits and bytes in different languages can be found in folder bits_manipulations.



