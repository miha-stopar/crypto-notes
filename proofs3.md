# Discrete logarithm equality 

We have Schnorr group G_q and elements g1, h1, g2, h2. Prover knows a secret x such that:

```
h1 = g1^x and h2 = g2^x
```

Chaum-Pedersen protocol [1] proves the knowledge of dlog and the fact that it is the same: log_g1(h1) = log_g2(h2).

![dlog equality](https://raw.github.com/miha-stopar/crypto-notes/master/img/dlog_equality.png)

# Blind non-interactive proof of discrete logarithm equality

Proof of dlog equality can be made non-interactive by using hash function - instead of challenge prepared by a verifier, the prover generates it as hash(A, B), where A = g1^r % p and B = g2^r % p are the two values sent in the first step to the verifier.

To make it blind,


# Pseudonym systems [2]

Lysyanskaya et al. presented in 1999 a pseudonym system [2] which allow users to interact with multiple organizations anonymously, using pseudonyms. There were similar systems (for example Chaum [3]) proposed before, but they did not protect the system against dishonest users who collectively use their pseudonyms and credentials (they share the identity). Also, previous systems relied heavily on a trusted center.

Chaum [3] introduced a system in 1985 unlinkable pseudonyms (or nyms) where organizations cannot extract any knowledge on the user even if they combine their databases - a user can obtain a credential from one organization using one nym and demonstrate possession of the credential to another organization with a second nym, where the second organization cannot extract the first nym even if it colludes with a first organization. For example, Bob may get a credential asserting his good health from his doctor (who knows him by one nym), and show this to his insurance company (who knows him by another nym).

To motivate the user not to share his pseudonyms or credentials with other users, Lysyanskaya et al. [2] proposed to have a certificate authority (CA) to enable a user to prove to an organization that his pseudonym actually corresponds to a master public key of a real user - if the user shares a credential to that pseudonym, his master secret key can be computed.

A similar system was proposed by Caneti et al. [5], however the difference (to [2]) is that here organizations grant credentials only to users who reveal their identity to them (and then these credentials can be used anonymously to interact with other organizations).

## Registration with CA

Each user must first register with the CA, revealing his true identity and his master public key, and demonstrating possession of the corresponding master secret key. A user's pseudonym with the CA is his master public key.

A user thus obtains a credential from CA which is used in the registration with an organization (next step). This credential tells that the user has a valid master secret and public key (certifying that he has a valid identity).

## Registration with an organization (pseudonym generation)

User U's master secret key is x and master public key is g^x.

![nym generation](https://raw.github.com/miha-stopar/crypto-notes/master/img/nym_generation.png)

At the end user U and organisation O remembers U's nym: (a, b). 

However, how is it ensured that the user doesn't use some random value for x (and not the secret that has been used for registration with CA)? Note that a is generated in the way that user does not know log_g(a). This is supposed to ensure that all user's nyms are tied to his secret x, but I think it does not.




# Non-transferable anonymous credentials [4]

While [2] discourage users from sharing their pseudonyms and credentials using PKI-assured non-transferability (sharing a credential implies also sharing a particular, valuable secret key from outside the system), Camenisch-Lysyanskaya [4] proposed a system which uses all-or-nothing non-transferability which does not require having some external valeable secret key. Here, sharing one pseudonym or credential implies sharing all of the user's other credentials and pseudonyms in the system.





[1] D. Chaum and T. P. Pedersen, Wallet databases with observers, Advances in Cryptology — CRYPTO ’92 (E. F. Brickell, ed.), LNCS, vol. 740, Springer-Verlag, 1993, pp. 89– 105.

[2] A. Lysyanskaya, R. Rivest, A. Sahai, and S. Wolf. Pseudonym systems. In H. Heys and C. Adams, editors, Proceedings of SAC 1999, volume 1758 of LNCS, pages 184–99. SpringerVerlag, Aug. 1999.

[3] David Chaum. Security without identification: transaction systems to make Big Brother obsolete. Communications of the ACM, 28(10), 1985.

[4] Camenisch, J. and Lysyanskaya, A., “Efficient non-transferable anonymous multi-show credential system with optional anonymity revocation,” in EUROCRYPT 2001, 93-118. 2001.

[5] Ran Canetti, Moses Charikar, Ravi Kumar, Sridhar Rajagopalan, Amit Sahai, and Andrew Tomkins. Non-transferable anonymous credentials. Manuscript, 1998. Revision in submission, 1999.


