These are by no means rigorous definitions, just brief sketches.

For notes below:
(E, D) is a cipher over (K, M, C), where E is encryption function, D decryption function, K key space, M message space, C cipher space, and D(k, E(k, m)) = m for all m from M and k from K.

# Perfect secrecy

Perfect secrecy is the strongest notion of cryptanalytic difficulty: for a given ciphertext c and two messages m_0, m_1 (the same length), the probabilities Pr[E(k, m_0)=c] and Pr[E(k, m_1)=c] are the same.

Basically, encrypted message provides no information about the original message.
That means no ciphertext-only attacks exist. 
For example One-Time Pad (OTP) has perfect secrecy. 
However, it can be shown that no cipher with keys shorter than plaintext can be perfectly secret - which is impractical because the two parties that communicate need to have a way to agree on really long keys (at least as long as the message) and such a mechanism could be then used to transmit the messages themself too.

# Pseudorandom generator (PRG)

PRG is deterministic procedure that maps a random seed to a longer pseudorandom string.
`G:{0, 1}**s -> {0, 1}**n, n >> s`

G is mapping seeds (`{0, 1}**s` is called seed space) which could be only 128 bits long to a much larger string which might be gigabytes long.

The idea is to replace random key by pseudorandom key - for example to make OTP practical (having longer pseudorandom key generated out of a shorter key).

PRG must be unpredictable - you cannot predict the next bit, even if you know all the previous ones.

# Stream cipher

Stream ciphers combine (xor) plaintext with pseudorandom cipher digit stream (keystream). Pseudorandom stream is generated using PRG.
This is the same approach as in OTP - but OTP uses a keystream of completely random digits, while stream ciphers use a small size key (like 128 bits) to generate a pseudorandom keystream (stream ciphers make OTP practical).

Here, security depends on the PRG that is used to generate pseudorandom keystream. And obviously - the same keystream must never ever be used twice (as in OTP), but keystream here is much longer than key and can be used to encrypt plaintext way longer that is the key length.

Stream cipher:

 * G is PRG
 * c = E(k, m) := m xor G(k)
 * D(k, c) := c xor G(k)

Stream cipher cannot have perfect secrecy because the key is (much) shorter than the message. A different definition of security is needed to argue that stream cipher is secure - semantic security.

# Semantic security

Semantic security is weaker than perfect secrecy. It implies that any information revealed from the ciphertext cannot be feasibly extracted. This is the computational analogue to perfect secrecy.
As in perfect secrecy definition, the ciphertext should not reveal any information about the original message (except the length), but reveal here means to be possible to write a computer program able to reveal the information.

The formal definition of semantic security uses a challenger / adversary game (such games are used in many security definitions):

 * adversary creates two messages of the same length and sends them to the challenger
 * challenger generates a random key, encrypts one of these two messages and sends the ciphertext to the adversary
 * adversary guesses which message was sent back

![semantic security](https://raw.github.com/miha-stopar/crypto-notes/master/img/semantic_security.png)

The accuracy of the guesses must not be better than a coin flip.

(One-time) semantic security guarantees secrecy for stream ciphers, where key is used only once (the adversary does not see two or more ciphertexts encrypted with the same key). For block ciphers security for many-time key needs to be defined - semantic security under the chosen-plaintext attack.

For example ECB mode is not semantically secure because it always maps identical plaintexts to identical ciphertexts. 
The attacker can create two messages, both of two blocks length. The first message should consist of two identical blocks, the second of two different blocks. The attacker can simply check whether the two blocks in returned ciphertext are the same or not.

# Secure Pseudo Random Functions (PRFs) and Pseudo Random Permutations (PRPs)

PRFs and PRPs are an abstraction of the concept of block ciphers.
The goal of the PRF is to look like a random function from X to Y.  

PRF is a function `F: K x X -> Y`, such that exists 'efficient' algorithm to evaluate `F(k, x)`.

PRP is a function `E: K x X -> X` such that:

 * exists 'efficient' deterministic algorithm to evaluate `E(k, x)`
 * the function (k, .) is one-to-one
 * exists 'efficient' inversion algorithm `D(k, x)`

Funs[X, Y]: the set of all functions from X to Y
S_F = {F(k, .) such that k from K} (subset of Funs[X, Y] - it is a set of all functions from X to Y defined by a given PRF F)

F is secure PRF if: 

 * the challenger chooses either k from K and have f = F(k, .) or truly random function f from Funs[X, Y]
 * the adversary sends x_1,...,x_q to the challenger
 * the challenger returns f(x_1),...,f(x_q)
 * the challenger cannot say from which set the function has been chosen

This means a random function from Funs[X, Y] is indistinguishable from a random function in S_F.

For secure PRP you replace the Funs[X, Y] with Perms[X].

For example 3DES and AES are believed to be secure PRPs (e.g AES-128: K x X -> X where K = X = {0, 1}^128)

The goal is to build a secure block cipher encryption from a secure PRP (like from AES). The problem is that here we have a many-time key (the same key for multiple blocks) - see semantic security under chosen-plaintext attack below.

# Chosen-plaintext attack (CPA)

This is an attack model which presumes that the attacker can choose arbitrary plaintexts and obtain the corresponding ciphertexts (and thus to gain some information which reduces the security of the encryption scheme). As a side note - an attack where the attacker knows a set of plaintexts and corresponding ciphertexts is called known-plaintext attack (the attacker does not have a power to choose which plaintext/ciphertext pairs are known).

An example might be a file storage system where the same AES key is used to encrypt files for many users (such a system would be terrible and should never be used). If Eve and Alice both use this system, Eve might start putting arbitrary files on the system and checks whether the resulting ciphertexts are the same as those of Alice (assuming that Eve can see Alice ciphertexts by eavesdropping or later by stealing the laptop). This way Eve could succeed in generating a plaintext which is encrypted in the same ciphertext as one of the Alice plaintexts.

# Semantic security under chosen-plaintext attack

A definition is similar to the (one-time) semantic security definition, but here an adversary sends `q` pairs of plaintexts (instead of one pair), the challenger encrypts either the first or the second messages and sends them to an adversary. Semantic security means an adversary cannot say which messages (first or second) were encrypted.

An adversary is much more powerful here, because he can use the knowledge gained from q_1,...,q_i to cleverly construct q_(i+1) and thus somehow win the game. Consider the following case.

An example of not CPA-secure system is a block encryption where CBC mode with predictable initialization vector is used. Here an adversary can win a CPA game in two steps:

 * first ask a challenger for an encryption of the next predicted initialization vector iv_1: E(k, 0 \xor iv_1)
 * let us denote the next predicted initialization vector as iv
 * send the two following messages to the challenger: m_0 = iv \xor iv_1 and an arbitrary m_1 which is different from m_0 
 * challenger returns either E(k, (iv \xor iv_1) \xor iv) or E(k, m_1). The first one is the same as E(k, iv_1) which is a known value from the first step, thus the adversary can easily say which message (m_0 or m_1) has been encrypted

![semantic security many time key](https://raw.github.com/miha-stopar/crypto-notes/master/img/semantic_security_many_time_key.png)

The famous BEAST attack exploits the misuse of initialization vectors as well, see more on http://crypto.stackexchange.com/questions/1078/how-can-cipher-block-chaining-cbc-in-ssl-be-attacked/1082#1082.

The adversary in CPA definition is so powerful that the encryption must be randomized to withstand the attack. If the same message is encrypted twice, the obtained ciphertexts need to be different (non deterministic encryption).

Two security notions have been summarized:

 * semantic security against one-time CPA
 * semantic security against many-time CPA

But neither ensures data integrity.

# Message integrity (MACs)

While integrity (and confidentiality) is a crucial part of authenticated encryption (see below the section about CCA-security), sometimes only integrity is needed, not confidentiality. For example to protect public binaries on disk or banner ads on web pages.

Alice sends message m to Bob. Alice needs to generate a tag S(k, m) which is later verified by Bob: V(k, m, tag) returns yes (the message m has not been changed) or no (the message m has been changed).

![mac](https://raw.github.com/miha-stopar/crypto-notes/master/img/mac.png)

The integrity is provided by MAC.
MAC I = (S, V) defined over (K, M, T) is a pair of algorithms:

  * S(k, m) outputs t in T
  * V(k, m, t) output 'yes' or 'no'

Integrity requires a secret key, otherwise attacker can easily modify message and re-compute the tag.

## What is a secure MAC?

The attacker can obtain tags for arbitrary messages of his choice (Alice sends m_1, m_2, ..., m_q to Bob and retrieves t_i = S(k, m_i) for each i). This is called a chosen message attack.

Having a secure MAC, the attacker should not be able to produce a valid tag for a new message (producing it is defined as existential forgery). 
Also, for a given (m, t) the attacker should not be able to produce (m, t1) for t != t1.

MAC is secure if it is existentially unforgeable under chosen message attack.

MAC I = (S, V) is secure (once again, this is not a formal definition and it omitts half of the stuff) when the probability that the adversary can get (m, t) such that `V(k, m, t) = 'yes'` and `(m, t) is not in {(m_1, t_1),...,(m_q, t_q)}` is negligible.

Any secure PRF gives a secure MAC:

For a PRF F: K x X -> Y we define a MAC I_F = (S, V) where:

 * S(k, m) := F(k, m)
 * V(k, m, t): output 'yes' if t = F(k, m) and 'no' otherwise

That means we compute tag simply as applying PRF on the message.

But this MAC is secure only if the |Y| is big enough - if for example Y = {0, 1}^10 that the probability that the adversary guesses the tag is not negligible.

Theorem: If F: K x X -> Y is a secure PRF and 1/|Y| is negligible then I_F is a secure MAC (like |Y| = 2^(80)).

So we can use AES for computing MACs for 16-byte messages, but how to convert small-MAC into a big-MAC (we do not want to have tags of the same length as messages when these are long)?

Two most used constructions to convert a small-PRF into a big-PRF are:

 * CBC-MAC (used in ANSI X9.9, X9.19, FIPS 186-3)
 * HMAC (used in SSL, IPsec, SSH)

# Collision resistance

Collision resistance plays an important role in providing message integrity. A very popular and widely used MAC called HMAC is built from collision resistant hash function.

Let H: M -> T be a hash function (|M| >> |T|). A collision for H is a pair m_0, m_1 such that:
 
`H(m_0) = H(m_1) and m_0 != m_1`

A function H is collision resistant if the probability that collision is found is negligible. A collision resistant function is for example SHA-256 (outputs 256 bits).

If I = (S, V) is a secure MAC for short messages (like AES) over (K, M, T) and H is collision resistant (H: M_big -> M), then I_big as defined below is secure:

I_big = (S_big, V_big) over (K, M_big, T) where:

```
S_big(k, m) = S(k, H(m))
V_big(k, m, t) = V(k, H(m), t)
```

This simply means you can use the hash value of the message to calculate a MAC and it will be a secure MAC if the hashing function is collision resistant.

But, can we use collision resistant hash function to directly (without PRF) build a MAC?

The first thing that might come to mind is to define a MAC as:

`S(k, m) = H(k || m)`

However, this is terribly insecure due to the inner construction of collision resistant hash functions (which all follow Merkle-Damgard paradigm, see about length-extension attack in [attacks.md](https://github.com/miha-stopar/crypto-notes/blob/master/attacks.md).

There is a standard method to convert a collision resistant hash function into a MAC, named HMAC:

`HMAC: S(k, m) = H(k xor opad || H((k xor ipad) || m))` 

where H is for example SHA-256; opad and ipad are fixed values.

## Generic attack on collision resistant functions

Let H: M -> {0, 1}^n be a hash function (|M| >> 2^n)

Generic algorithm to find a collision in time O(2^(n/2)):

 * Choose 2^(n/2) random messages in M
 * For i = 1,...,2^(n/2) compute t_i = H(m_i)
 * Look for a collision (t_i = t_j). If not found, go to the first step.

The expected number of iterations is about 2.

This on the first glance somehow surprising result is due to the so-called birthday paradox: If r_1,...,r_n are randomly chosen from {1,...,B}, then when n = 1.2 x B^(1/2) it holds: Pr[exists i != j: r_i = r_j] >= 0.5.

function | digest size (bits) | speed (MB/sec) | generic attack time
-------- | ------------------ | -------------- | -------------------
SHA-1    | 160                | 153            | 2^(80)
SHA-256  | 256                | 111            | 2^(128)
SHA-512  | 512                | 99             | 2^(256)
Whirlpool| 512                | 57             | 2^(256)

# Chosen-ciphertext attack

Attack model which presumes that the attacker can choose a ciphertext and obtain its decryption (by somehow tricking the system) under an unknown key (and thus potentially recover the key) or at least some partial information about the plaintext (for example an attacker can sends TCP/IP packet and retrieves ACK from the server, in this case the attacker knows that the checksum is valid).

Note that in a famous CBC padding oracle attack (Matasano challenge 17) there is not even plaintext needed - you just feed in the ciphertext and use the information whether it can be decrypted to something with valid padding. However, the key is not retrieved, only the plaintext.

A beautiful chosen-ciphertext attack is Bleichenbacher attack against PKCS#1 (Matasano challenge 48).

To summarize: semantic security which was defined by Goldwaser and Micali [1] should guarantee that an adversary is not able to obtain any partial information about a message given its encryption.
However, this guarantees security only when the adversary is passive.

If adversary is able to mount an active attack, we need a different notion of security. Rackoff and Simon [2] defined security against adaptive chosen ciphertext attack - if an adversary can inject messages into a network (may be encryptions), the adversary might be able to extract partial information about the cleartexts.

Rackoff and Simon allowed an adversary to obtain decryptions of his choice (decryption oracle), but the adversary still should not be able to extract any partial information about the cleartexts (adversary is not allowed to pass to the oracle the ciphertext of the message that he is trying to decrypt, but he can pass any other). See next section for some details.

A nice overview of the security definitions is given in Cramer and Shoup [3].

# Semantic security under the chosen ciphertext attack (CCA)

The CCA security is required against an active attacker who can tamper the messages. For example CBC mode (even with random initialization vector) is not semantically secure under CCA.

CPA-secure system is secure against eavesdropping only. Integrity is needed as well to prevent active attacks.

## Ciphertext integrity

Ciphertext integrity game goes roughly as follows:

 * challenger chooses a random key k
 * an adversary sends m_1,...,m_q messages to the challenger
 * the challenger returns encrypted messages to the adversary: c_1 = E(k, m_1),...,c_q = E(k, m_q)
 * the challenger creates a new ciphertext c (which is not from {c_1,...,c_q})

(E, D) has ciphertext integrity if the probability that an adversary can create a ciphertext that is not rejected is negligible.

## Chosen-ciphertext security

Adversary have both CPA and CCA power:

 * can obtain the encryption of arbitrary messages of his choice
 * can decrypt any ciphertext of his choice, other than challenge
 * challenger chooses a random key k
 * adversary submitts queries to the challenger - every query can be either CPA query (submitts two messages and receives the encryption of either first or second - these ciphertexts c_i are called challenged ciphertexts) or CCA query (submitts an arbitrary ciphertext and gets back its decryption); but the ciphertext submitted in CCA query is not allowed to be one of the ciphertexts obtained by CCA query (otherwise the adversary would easily win) 

Note that CPA and CCA queries can be interleaved in any order. When asking for a decryption for c_i (CCA query), c_i must not be from {c_1, ..., c_(i-1)} where c_1,...,c_(i-1) are all responses to CPA queries so far.

![semantic security](https://raw.github.com/miha-stopar/crypto-notes/master/img/semantic_security.png)

The system is chosen-ciphertext secure (as always here - this is just a sketch of definition) if the probability that the adversary determines which messages have been returned in CPA queries is negligible (the attacker is allowed to decrypt any ciphertext other than the challenged ciphertexts, but he still cannot determine whether the first or second messages have been encrypted).

(for pub-key is different, see [public_key_security.md](https://github.com/miha-stopar/crypto-notes/blob/master/public_key_security.md))
Note that for public key encryption the definition for security is different (no need to give an attacker the ability to mount the CPA, because an attacker has a public key and can encrypt any message), see [public_key_security.md](https://github.com/miha-stopar/crypto-notes/blob/master/public_key_security.md).

## Authenticated encryption

Authenticated encryption (E, D) is a cipher where:

 * E: K x M (x N) -> C (N is a nonce space)
 * D: K x C (x N) -> M U {⊥}

Decryption returns either decrypted message or bottom symbol (⊥) when cipher text is rejected.

Authenticated encryption must provide:

 * semantic security under a CPA
 * ciphertext integrity (attacker cannot create new ciphertexts that decrypt properly)

## Conclusion

Why is authenticated encryption so important? Because it can be proved that authenticated encryption provides CCA security. Thus authenticated encryption ensures confidentiality against an active adversary that can decrypt some ciphertexts (but it does not prevent replay attacks or side channels attacks).

# Authenticated encryption constructions

As there are CPA-secure systems and secure MACs - how can be an authenticated encryption build out of this? 
Authenticated encryption has been introduced in 2000, before this developers had to combine API for CPA-secure encryption (like CBC with random IV) and API for MAC (like HMAC). However, many mistakes can be made when combining these two APIs (the proper way is encrypt-then-MAC).

Today there are modes of operation which provide authenticated encryption, like:

 * GCM
 * CCM
 * EAX

All these three modes are nonce-based (nonce has to be unique per key, but it does not to be random). Also, all these three modes support authenticated encryption with associated data (AEAD) - not the full message is intended to be encrypted, but the full message has to be authenticated (like network packet - header which contains the destination is not encrypted).

# Standard model and random oracle model

When some scheme or protocol is proved to be secure, it means that it is secure under some assumptions. In standard model (also plain or bare model) the assumption is that the adversary is only limited by the amount of time and computational power.

As achieving security proofs in a standard model are difficult, in many proofs cryptographic primitives are replaced by idealized versions. In random oracle model a hash function is replaced with a random function.

However, a truly random function is slow and we would need to allow all parties involved in the security proof to be able to run exponential number of steps. This means also adversary would not be given only a limited computational power. To avoid this, hash operations in the scheme/protocol are replaced with an external call to a random oracle to which all parties have access (and oracle calls are assumed to take no time). Thus, the parties are still limited to a polynomial time computations.

However, hash functions are not random functions and it has been shown that there are protocols which are secure in the random oracle model but insecure in standard model.




[1] S. Goldwasser and S. Micali. Probabilistic encryption. Journal of Computer
and System Sciences, 28:270–299, 1984.

[2]  C. Rackoff and D. Simon. Noninteractive zero-knowledge proof of knowledge
and chosen ciphertext attack. In Advances in Cryptology–Crypto
’91, pages 433–444, 1991.

[3]  R. Cramer and V. Shoup, A practical public key cryptosystem provably secure against adaptive chosen ciphertext attack, Advances in Cryptology – Crypto ’98, LNCS 1462, Springer-Verlag, 1998, pp. 13-25.



