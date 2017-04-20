# Commitment schemes

Let's denote a commitment to a value v as C = Comm(v). The recipient cannot infer any information about v given C (property called hiding) and cannot find v1 which is different from v and Comm(v1) = C (this means the prover cannot change her mind later, the property is called binding).

Commitment can be thought of as an envelope - the value is put in the envelope and sealed. Upon receiving a commitment, a verifier cannot tell what is inside (hiding) and the prover cannot sneak another value in the envelope later (binding).

A commitment can be later open to reveal what is inside.

Commitment schemes have many applications. Consider for example Alice and Bob flipping the coin and communicating over the phone (they are not in the same place). Instead of reporting the actual coin flip, Alice only tells Bob a commitment of the value. Bob then flips the coin and reports the result. Alice now reveals what she commited to (Bob cannot cheat).

For commitments schemes the following two properties need to hold (informally):

 * secrecy (hiding): (almost) no information is revealed to the verifier
 * unambiguity (binding): your chances of being able to change your mind are very small (the sender cannot change the message)

There is a special subset of commitment schemes called trapdoor commitment schemes. Trapdoor commitment schemes allows to overcome the binding property if you know a special information (called trapdoor).
This means that a committer can change her mind and open a commitment abiguously.

**We will see below how trapdoor commitment schemes can be used to build zero knowledge proofs out of sigma protocols. Note: it is much easier to prove that a protocol is a sigma protocol than to prove that a protocol is a zero knowledge proof of knowledge.**

## Pedersen commitment scheme

Receiver chooses:

 * large primes p and q such that q | p-1.
 * generator g of order q subgroup of Z_p*
 * random secret a from Z_q
 * h = g^a % p

Values p, q, g, h are public, a is secret (a is called trapdoor).

When committing to some x from Z_q, sender chooses random r from Z_q and sends `c = (g^x * h^r) % p` to receiver.

![pedersen commmit](https://raw.github.com/miha-stopar/crypto-notes/master/img/pedersen_commit.png)

When a sender reveals x and r, receiver verifies that `c = (g^x * h^r) % p`.

![pedersen decommit](https://raw.github.com/miha-stopar/crypto-notes/master/img/pedersen_decommit.png)

Note that Pedersen commitment is a **trapdoor commitment scheme**: if the committer knows the exponent a (trapdoor), she can decommit c to an arbitrary value x'. 
Decommitter simply reveals (x', r' = r + (x - x')/a). It holds:

```
(g^x' * h^r') % p = (g^x' * h^r * h^((x-x')/a)) % p = (g^x' * h^r * g^(x-x')) % p = g^x * h^r = c 
```

# Zero-knowledge proofs

The majority of the notes in this section are taken from Hazay-Lindell [1] Sigma protocols and efficient zero-knowledge chapter which is based on [2].

Zero-knowledge proof is an interactive proof where verifier learns nothing beyond the correctness of the statement being proved. Zero-knowledge proofs can be used to enforce the proper behavior of the verifier.

Many zero-knowledge protocols are constructed from a simpler primitive called sigma protocol (it is far easier to construct a protocol and prove that it is a sigma protocol, than to construct a protocol and prove that it is zero-knowledge proof of knowledge).

## Relation to interactive proofs

In the interactive proof the prover wants to ascertain a verifier whether a given string belongs or not to a formal language.
The prover posseses unlimited computational resources, but cannot be trusted. The verifier has limited computation power.

For the interactive proof two requirements need to hold:

 * completeness: if the string belongs to a language, the honest verifier (the one that follows the protocol properly) can be convinced of this fact by un untrusted prover
 * soundness: if the string does not belong to the language, no prover can convince the honest verifier that it belongs, except with some small probability

Zero knowledge proof is an interactive proof with the additional property that the verifier learns nothing but the correctness of the statement (that the string belongs to the language).

## Schnorr's protocol

A simple example of a zero knowledge proof is the proof of knowledge of a discrete logarithm. This can be done by Schnorr's protocol:

 * prover wants to prove that t = g^s % p
 * system parameters are large primes p and q such that q | p-1; g is a generator of a subgroup of order q in Z_p* (subgroup is Schnorr group, see [number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md))
 * prover chooses random r from Z_q and sends x = g^r % p to verifier
 * verifier knows t and wants a proof that prover knows its discrete logarithm (**t can be thought as a public key, while s is a secret key**)
 * verifier chooses random c from [1,...,2^n] and sends it to a prover (n is fixed such that 2^n < q)
 * prover sends y = (r+s*c) % q to the verifier
 * verifier verifies that g^y = g^r * (g^s)^c (mod p)

![schnorr protocol](https://raw.github.com/miha-stopar/crypto-notes/master/img/schnorr_protocol.png)

Note that a verifier cannot extract s from y = r + s*c since she does not know r (to extract r a discrete logarithm would need to be solved). Also note that if c is not chosen randomly (if prover can predict the value c), the prover can set x = g^y * t^(-c) and convince the verifier that she knows the secret witout really knowing it.

At the end of the Schnorr's identification protocol the verifier is convinced he is interacting with the prover, or better to say - with someone who knows the secret key that corresponds to the prover's public key. Which means with someone who knows the dlog of t.

Schnorr's protocol can be used as identification protocol and is different from the standard password authentication where prover's secret key is her password pw and public key is hash(pw). Here the prover sends hash(pw) and the verifier checks whether this matches the stored value. However this is not secure agains eavesdropping attacks as the adversary who observes one interaction between prover and verifier and sees hash(pw), can later impersonate a prover.

### Is Schnorr's protocol a zero knowledge proof?

The required properties for zero knowledge proofs are:

 * completeness: if the statement is true, the honest verifier (the verifier that follows the protocol properly) will be convinced of this fact with overwhelming probability
 * soundness: no one who does not know the secret can convince the verifier with non-negligible probability
 * zero knowledge: the proof does not leak any information

Zero knowledge property is usually demonstrated with an existance of a simulator that produces a fake transcript of protocol messages that is indistinguishable from the actual protocol messages. If all messages can be thus simulated, verifier does not learn anything new during the protocol.

There is also a notion of honest-verifier zero knowledge protocol which means that the protocol is zero knowledge if the verifier follows the protocol honestly.

**Why is the size of challenge (t) important - because if the prover can guess the challenge, she chooses z and sets the first message as g^z * h^(-e) - and thus convince the verifier.**

Soundness can be demonstrated by an existance of a knowledge extractor algorithm which can extract the secret from a prover. The intuition behind is - if a prover can convince the verifier that she knows the secret without actually knowing it, then no algorithm could extract the secret from a prover.

The Schnorr's protocol is sound because there exists an algorithm which extracs a secret from a prover. If prover sends x to the verifier and then follows the protocol two times (obtains two challenges c_1 and c_2) and returns proper y for both cases (y_1 = r + s\*c_1 and y_2 = r + s\*c_2), this means that x = g^(y_1)*t^(-c_1) = g^(y_2)*t^(-c2) from where it follows s = (y_2 - y_1)/(c_2 - c_1) - prover knows a secret s. This is called rewinding.

Schnorr's protocol is however honest-verifier zero knowledge (below more about this), it is not known if it is zero knowledge.

Zero-knowledge property can be shown by the existance of a simulator with which we can simulate acceptable protocol runs without knowing the actual language statement.

This can be done for Schnorr protocol if t is small (there are no many possible challenges c). Here again, a rewinding technique is used. We choose y1, try to guess c1, compute x1 = g^y1 * t^(-c1). If we now send x1 to the verifier and get c1 (for which there is a high probability, because the space for challenges is small), the transcript (x1, c1, y1) will be accepted. So, we can simulate acceptable runs without knowing s (s is secret, we know only t = g^s % p).

As mentioned for large t, we don't know whether Schnorr protocol is zero-knowledge (whether verifier could extract some information about s after having a number of transcripts). However, a description how to make sigma protocols zero-knowledge is given below.

## Schnorr signatures

Schnorr signature is obtained by applying the Fiat-Shamir transformation to the Schnorr protocol.

The setup of parameters is as in Schnorr id scheme described above. The difference is that instead of challenge c being chosen by a verifier, prover calculates c = hash(m || g^r), where m is a message that is to be signed.
The signature is then (y, c), where y = r + s\*c

Verifier needs to check g^y = g^r * (g^s)^hash(m || x) (mod p).

![schnorr signature](https://raw.github.com/miha-stopar/crypto-notes/master/img/schnorr_signature.png)

## Formal definition of sigma protocols

Let R be a binary relation, where the only restriction is that if (x, w) is from R, then the length of w is at most p(|x|) for some polynomial p(). We may think of x as an instance of some computational problem.
For example for discrete log problem we can define relation: 

```
R_DL = {((G, q, g, h), w); h = g^w}
```

where G is group of order q with generater g, q is prime, and w is from Z_q.
R_DL contains the entire set of discrete log problems with their solutions.

We consider the protocols of the following form:

 * Common input: the prover P and verifier V both have x.
 * Private input: P has a value w such that (x, w) is from R.
 * The protocol template:
	* P sends V a message a.
	* V sends P a random t-bit string e.
	* P sends a reply z, and V decides to accept or reject based solely on the data it has seen (based on value x, a, e, z).

A transcript (a, e, z) is an accepting transcript for x if the protocol instructs V to accept based on the value (x, a, e, z).

A protocol is sigma protocol for relation R if it is a three-round public-coin protocol of the form above and the following requirements hold:

 * Completeness: If P and V follow the protocol on input x and private input w to P where (x, w) from R, then V always accepts.
 * Special soundness: There exists a polynomial-time algorithm A that given any x and any pair of accepting transcripts (a, e, z), (a, e', z') for x, where e!=e', outputs w such that (x, w) from R.
 * Special honest verifier zero knowledge: there exists a probabilistic polynomial-time simulator M, which on input x and e outputs a transcript of the form (a, e, z) with the same probability distribution as transcripts between the honest P and V on common input x.

It is clear that Schnorr's protocol is sigma protocol - the simulator for special soundness is described above (you have two transcripts (a, e, z) and (a, e', z'), so you have g^z = g^r * g^(e * w) and g^z' = g^r * g^(e' * w), from where you get g^((z-z')/(e-e')) = g^w and you have w). The simulator for special honest verifier zero knowledge is: you take random z and random e, compute a = g^z * h^(-e) (you don't need to know w) and you have a transcript which has exactly the same probability distribution as a real conversation between the honest prover and honest verifier.
However there are many other examples of sigma protocols.

Note that sigma protocol also proves that x has a solution (in dlog there is always a solution, but in other cases this is not necessary) with probability > 1-2^(-t).
This is because if the solution w does not exist, there exists for every a at most one query e that results in an accepting transcript (otherwise w would be obtained due to special soundness property and thus the solution of x would exist), so the prover would convince the verifier with the probability at most 2^(-t).

## Sigma protocol for Diffie-Hellman tuples

Given a tuple (G, q, g, h, u, v) can we prove that there exists w such that u = g^w and v = h^w?
Let's define a protocol:

 * P and V both have (G, q, g, h, u, v)
 * P has a value w such that u = g^w and v = h^w

The protocol:

 * P chooses a random r from Z_q, computs a = g^r and b = h^r and sends (a, b) to V.
 * V chooses a random challenge e from {0,1}^t where 2^t < q and sends it to P.
 * P sends z = r + e * w mod q to V, which checks that g^z = a * u^e and h^z = b * v^e, that G is a group of order q wither generators g and h, and that u, v are from G.

If we show that this is sigma protocol, we have a way to prove that there exists w such that u = g^w and v = h^w.
The completeness holds because (similarly for h^z = b * v^e):

```
g^z = g^(r + e * w) = g^r * g^(w * e) = g^r * (g^w)^e = a * u^e
```

Special soundness: for two transcripts ((a, b), e, z) and ((a, b), e', z') we get w as for the Schnorr's protocol from either of the two equations: g^(z-z') = u^(e-e') and h^(z-z') = v^(e-e').

The special honest verifier zero knowledge property can be shown as for the Schnorr's protocol.

## Some properties of sigma protocols

**If you execute sigma protocol in paralell l-times, the resulting protocol is a new sigma protocol with challenge length l * t.**

Note that the probability that a prover guesses a challenge is 2^(-t) for one execution, if there are l executions, the probability is 2^(-l * t). Thus, the soundness error is reduced to 2^(-l * t).

**If there exists a sigma protocol SP1 for R with challenge length t, then there exists a sigma protocol SP2 for R with challenge length t', for any t'.**

If t' is smaller than t, we define a new sigma protocol SP2 as follows - SP2 sends the first message as in SP1. V sends a random t'-bit string e. In the last step SP2 appends t-t' 0s to e and use this value to compute the message z. V then checks z as in SP1.
Completeness and special soundness obviously still holds. Special zero knowledge holds because the simulator M (which exists for SP1) works for every challenge e (it generates an accepting transcript out of e).

If t' is bigger than t, the SP2 can be constructed by running SP1 in parallel multiple times and then adjust the length down (as above) to t' if needed.

## Proofs of knowledge

We haven't proved that sigma protocol is really a proof of knowledge yet. 
Intuitively, the special soundness is equivalent to the proof of knowledge, however, the proof should say how the two accepting transcripts are found.

The standard definition of proofs of knowledge which is given below was proposed by Bellare and Goldreich [4]. Roughly, some older definitions went like this: if the verifier is always convinced, then there exists a knowledge extractor which can extract a knowledge from a prover.
Bellare and Goldreich proposed definition which takes into account cases when the prover is able to convince the verifier with some non-constant probability smaller than 1.

Let's first formally define what is a proof of knowledge.
Let kappa: {0, 1}* -> [0, 1] be a function. A protocol (P, V) is a proof of knowledge for the relation R with knowledge error kappa if the following two properties are satisfied:

 * Completeness: If P and V follow the protocol on input x and private w to P where (x, w) from R, then V always accepts.

 * Knowledge soundness (or validity): There exists a constant c > 0 and a probabilistic oracle machine K (knowlege extractor) such that for every interactive prover function P\* and every x from L_R (L_R are all x for which a solution w exists), the machine K satisfies the following condition.
Let epsilon(x) be the probability that V accepts on input x after interacting with P\*.
If epsilon(x) > kappa(x), then upon input x and oracle access to P\*, the machine K outputs w such that (x, w) is from R within an expected number of steps bounded by:

```
|x|^c/(epsilon(x) - kappa(x))
```

**If SP is a sigma protocol for a relation R with challenge length t, then SP is a proof of knowledge with knowledge error 2^(-t).**

Completeness holds by definition. What we need to show is that we are able to find two accepting transcripts "sufficently" quickly.
Let's visualize all possible interactions in a matrix - one row for each possible choices by P, one column for each possible challenge value e.

![zero knowledge](https://raw.github.com/miha-stopar/crypto-notes/master/img/zero_knowledge.png)

I assume by this is meant that each row presents a possible first message in a protocol. 
The matrix contains 1 if P can convince V for the first message (which defines this row) and for the challenge that defines this column (there are 2^t possible values for challenge), otherwise there is 0.
The goal is to find two 1's in one row - when this is done, the special soundness gives as the result.

epsilon(x) is the probability that V is convinced (for each x there is its own matrix)- epsilon(x) is the fraction of 1's in the matrix, but this does not say anything about how 1's are distributed in rows.
If we start searching in a row where there is only a single 1, we have to break at some point and start looking in other rows (but we don't know if we are in a row like this).

It turns out we can design an algorithm which will find two 1's in a row "sufficiently" quickly.

## Proving compound statements

It is obvious how to construct a sigma protocol to prove an AND statement. Simply, the prover proves both statements in parallel - the verifier sends a single challenge for both proofs.

For OR statements the construction is more complex. Both, P and V, have a pair (x_0, x_1). P wants to prove to V that it knows w such that either (x_0, w) from R_0 or (x_1, w) from R_1 (w is solution for either x_0 or x_1).

Let b be such that (x_b, w) is from R_b. The protocol goes like this:

 * P computes the first message a_b using (x_b, w) as input. P chooses e_(1-b) at random and runs the simulator M on input (x_(1-b), e_(1-b)) - let (a_(1-b), e_(1-b), z_(1-b) be the output of M. P sends (a_0, a_1) to V.

 * V chooses a random t-bit string s and sends it to P.
 * P sets e_b = s xor e_(1-b) and computes the answer z_b to challenge e_b (the other z is from the output of M). P sends (e_0, z_0, e_1, z_1) to V.
 * V checks that e_0 xor e_1 = sand that both transcripts are accepting.

## Zero knowledge from sigma protocols

As we saw above sigma protocol is a proof of knowledge. However, it is not a zero knolwedge proof of knowledge - because it requires an honest verifier.

**The problem is that the challenge generation of the verifier cannot be simulated. 
Note that the point of zero knowledge proof is that the verifier learns nothing about the statement (except that it is correct), so if the verifier cheats and doesn't follow the protocol properly, it might be able to somehow extract some knowledge. 
For this reason we have to be somehow assured that verifier does not cheat.**

This can be solved by having the verifier commit to its challenge e before the execution begins (this way the verifier cannot choose some special e once he knows the first message a). The protocol then goes:

 * V chooses a random t-bit string e and interacts with P via the commitment protocol in order to commit to e.
 * P computes the first message a and sends it to V.
 * V decommits to e to P.
 * P verifies the decommitment and aborts if it is not valid. Otherwise, it computes the answer z and sends it to V.
 * V accepts if and only if transcript is accepting.

**If a commitment scheme is perfectly-hiding, then the sigma protocol extended as above is a zero knowledge proof for L_R with soundness error 2^(-t).**

The complexity of the protocol above depends on the cost of running a perfectly-hiding commitment protocol.
One possibility is to use Pedersen commitment protocol (see Commitment schemes section above) which is highly efficient.

However, the zero knowledge protocol above is not a proof of knowledge. The problem is that extractor cannot send prover two different e != e' for the same a, because it is committed to e which was sent in the first step.

Let's extend the protocol in the following way (we need a perfectly-hiding trapdoor commitment scheme):

* V chooses a random t-bit string e and interacts with P via the commitment protocol in order to commit to e.
 * P computes the first message a and sends it to V.
 * V decommits to e to P.
 * P verifies the decommitment and aborts if it is not valid. If valid, it computes the answer z. It sends z and trapdoor to V (verifier cannot use trapdoor to cheat anymore, because he already decommitted) .
 * V accepts if and only if the trapdoor is valid and the transcript is accepting.

The difference is that now the extractor can obtain the trapdoor and get e' != e (using trapdoor) for the same a. 

**It can be proved that this protocol is zero-knowledge proof of knowledge.**

![zero knowledge from sigma](https://raw.github.com/miha-stopar/crypto-notes/master/img/zk_from_sigma_protocol.png)

Note that Pedersen commitment is a trapdoor commitment scheme (see the first section on this page) and can be thus used to construct a zero-knowledge proof of knowledge.


[1] C. Hazay and Y. Lindell. Efficient Secure Two-Party Computation: Techniques and Constructions. Springer, 2010.

[2] I. Damgard. On Σ Protocols. http://www.cs.au.dk/~ivan/Sigma.pdf

[3] https://tools.ietf.org/html/rfc4252#page-8

[4] M. Bellare and O. Goldreich: On defining proofs of knowledge: proc. of Crypto 92.


