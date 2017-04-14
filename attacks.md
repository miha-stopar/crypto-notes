# Length-extension attack (Matasano challenge 29)

There is a common question (see security_definitions) - how to use hash function (usually function that does not require a key) to construct a secure MAC (keyed function).
It is a common mistake to use it this way:
insecureMAC(k, m) = H(k || m), where H is a hash function like SHA-256. 
This is terribly insecure. HMAC should be used instead. The problem is that the adversary can construct a new message and its hash without knowing the key. However, the adversary needs to know the key length.

Once the adversary knows the key length (he can guess it), he can work out the padding of the (key+message) which is used inside SHA1. He should just follow the SHA1 specifications (messag followed by 0s followed by a 64-bit integer to have a padded message of length 512 * n. The 64-bit integer should be the length of the original message - this is why you need the key length - as the original message is here key+message).

Now that an adversary has m + padding, he can add an arbitrary appendix (for example ";admin=true") and knowing the SHA1 internals he can compute the MAC of this constructed message (m + padding + appendix).

SHA1 uses Merkle-Damgard paradigm, meaning that the message is split into blocks and then the output of the previous block is passed as input into the compression of the current block. 

Value H_k is given as insecureMAC(k, m): 
H_k = h(H_(k-1), m[k-1]) = insecureMAC(k, m)

Value H_k can be then pushed as input to the compression function to calculate insecureMAC(message + padding + appendix):

insecureMAC(message + padding + appendix) = H_(k+1) = h(H_k, m[k] || PB) = h(insecureMAC(k, m), m[k] || PB) where PB is once again padding, but this time for the constructed message (note that the length here means (key+message+appendix).length).

Note: only the key length is required to be known.

# E = 3 RSA broadcast attack (Matasano challenge 40)

When the same message is sent to three different entities, each having its own RSA modulus N (let us say N1, N2, N3) and each using the same small encryption exponent e (like e = 3), the message can be trivially decrypted by an eavesdropper.

The message m is encrypted as `c1 = m^3 % N1, c2 = m^3 % N2, c3 = m^3 % N3`. The integers N1, N2, N3 are pairwise relatively prime (otherwise the factorization of N gets much easier, see RSA problems section in [public_key_security.md](https://github.com/miha-stopar/crypto-notes/blob/master/public_key_security.md)), and thus by Chinese remainder theorem we can find x such that:

 * x = c1 mod N1
 * x = c2 mod N2
 * x = c3 mod N3

Since m^3 < N1 * N2 * N3, x must be m^3 (Chinese remainder theorem says also that all solutions are congruent modulo the product N1 * N2 * N3).

Let us define N = N1 * N2 * N3. Using the Extended Euclidean algorithm we can find `r_i, s_i` such that `r_i * N_i + s_i * N / N_i = 1` (N and N/N_i are relatively prime).

 * `e_1 = s_1 * N_2 * N_3`
 * `e_2 = s_2 * N_1 * N_3`
 * `e_3 = s_3 * N_1 * N_2`

It holds:

 * `e_i = 1 mod N_i` (since `r_i * N_i + e_i = 1`)
 * `e_i = 0 mod N_j`, when j is not i

The solution x can be constructed as: `x = c1 * e_1 + c2 * e_2 + c3 * e_3`.

# Bleichenbacher's RSA signature forgery (Matasano challenge 42)

It works when exponent is 3. An RSA signature is pretty much the same as RSA encryption, but in a reversed process: the message is signed using the private key exponent and verified by the public key exponent.

First the message needs to be hashed. Let's say that we want to sign a message 'hi mom'. Its hash is: '925a89b43f3caff507db0a86d20a2428007f10b6'. The hash is then prepended with bytes in ASN.1 format which say what hash algorithm is used. We will simplify this and say only the hash length is prepended which is 20 (14 is hex). This is then prepended with PKCS-1 padding like this (byte of 0, byte of 1, string of 0xFF bytes, byte of 0):

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/bleichenbacher1.png

Note that n has to be big enough - in the case above the message has been formatted to have 384x8 bits = 3072 bits), so we have to choose p and q of 1536 bits to have n of 3072 bits.

On the message in this format is then applied the modular exponentiation (s = m ^ d % n) with a private key d. This is a signature of the message. The verifier needs to apply the modular exponentiation with a public key e. This means the verifier calculates s ^ 3 % n, checks the padding to extract the hash data, and finally verifies that the hash data is really the hash value of the message.

The problem appears when the implementation of the verifier does not check whether there is nothing behind the hash data. If the implementation checks for example only for a pattern like 'some string of 0xFF bytes and then a byte of 0' and then reads a specified number of bytes (given immediately after the pattern) to extract the hash data, the attacker can append a cleverly calculated string of bytes behind the hash data and provide a valid signature for it.

For example an attacker can provide a valid signature for 'hi mom' without knowing the private key exponent. First, the hash is padded like below:

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/bleichenbacher2.png

Only three 0xFF bytes are given before the byte of 0, byte 14 (hash length) and hash data. We assume here that implementation checks for pattern 'FF FF FF 00', then reads the hash length and finally extracts the hash data. The bytes behind the hash (garbage data) are never read (there are just bytes of 0 there, but this is to be changed immediately).

Let's denote the constructed string of bytes as x. Now, if the attacker defines garbage data of x so that it is a perfect cube (y ^ 3 = x), these means the valid signature has just been constructed (y).

Let's denote the constructed string with only bytes of 0 appended (as on the picture above) as x0. A cube root of x0 is then calculated or better to say its approximation is calculated (this can be done quickly using some nth root algorithm). Let's denote this value (cube root of x0) as r0. Now we have to check whether r0 ^ 3 > x0 (as said r0 is just an approximation of cube root). If it is not, we choose (r0+1) or something like this. We denote the selected value as r.

However, r ^ 3 should not be greater than the number representation of the string on the picture below (otherwise the bytes other than in garbage data will be changed): 

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/bleichenbacher3.png

The value r (or bytes representation of it) is the valid signature of 'hi mom' and can be sent to the verifier.

# Baby Bleichenbacher PKCS 1.5 Padding Oracle (Matasano challenges 46)

Let us assume that there is an oracle function on the server which is receiving RSA-encrypted messages and somehow returns whether the plaintext of the message is even or odd.

So, we have an encrypted message and no keys. We send an encrypted message to the server and we receive info about whether the plaintext (corresponding to the encrypted message) is even or odd. It turns out this is enough to decrypt the message. You should be aware that a plaintext is simply (transformed into) a huge number.

We use the oracle function iteratively and by each step we will have a more precise estimation in which range of numbers our plaintext (actually its corresponding number) is.

For an easier understanding below is an example with plaintext number being 7 and RSA n parameter being 11.

.. image:: https://raw.github.com/miha-stopar/crypto-notes/master/img/rsa_parity.png

We do not know the number 7, so let us say it is variable x. We can check whether 2*x is odd or even (actually whether (2*x)\*\*e is odd or even, but we can simply check 2\*x because (2\*x)^e is odd if and only if 2*x is odd). We know that n is odd, because it is a product of two odd (large prime) numbers. From this it follows that 2*x is odd when 2*x mod n is even. 

Thus, if we know that 2*x % 11 is odd, then x is bigger then 11/2. If 2*x % 11 is even, we know modulo 11 has not been applied and 2*x < 11.

We check 4*x then. If 4*x % 11 is even, it means it is one of the two following options:

 * 4*x < 11
 * 4*x >= 2 * 11 and 4*x < 3 * 11

In the first step (2*x % 11 is odd, so x >= 6) we saw x >= 6. So it is the second option (4*x >= 22 and 4*x <33). This means x >= 6 and x < 8.

We check 8*x then. We know 8*x % 11 is odd (56 % 11 = 1), so it is one of the following options:

 * 8*x > 11 and 8*x < 22
 * 8*x >= 33 and 8*x < 44
 * 8*x >= 55 and 8*x < 66
 
We know from the previos step that x is 6, 7, or 8. We know that 8*6 = 48 and this is not one of the three options. Thus, x is 7 or 8.

We repeat the process once more and we get: 7. It holds that when (2**i) * x is odd, we can estimate the minimum value, when it is even, we can estimate the maximum value. The code can be seen below:

```
  i = 1
  max = n
  min = 0
  while true do
    c1 = ((2**i)**e * c) % n
    isOdd = oracle.isOdd(c1)
    if isOdd
      min = max - (max-min)/2
      if min != max
        if ((min * (2**i)) % n) % 2 != 1
          min = min + 1
        end
      end
    else
      max = min + (max-min)/2
      if max != min
        if ((max * (2**i)) % n) % 2 == 1
          max = max - 1
        end
      end
    end
    if min == max
      break
    end
    i += 1
  end
  puts min # that is the value that is being sought
```

# Bleichenbacher PKCS 1.5 Padding Oracle (Matasano challenges 47, 48)

The Bleichenbacher attack is possible when the attacker has an access to the oracle that decrypts a ciphertext c into plaintext m and reveals whether m is PKCS compliant or not. 
That means whether: 

 * it starts with two bytes 00||02
 * the length of padding string (PS) is larger than 8
 * PS contains no zero bytes 00
 * PS is followed by a zero byte 00

Let us denote (N, e) the RSA public key and d as the corresponding secret key.
Encryption goes like:

```
c = m^e % N
```

Decryption:

```
m = c^d % N
```

Message m is compliant which means that it lies in the interval [2 * B, 3 * B - 1], where B = 2^(8 * (l-2)) and l is length of N in bytes (B is 00||01||00||...||00).

The attacker now submits many different values (starts with s = N/(3 * B))

```
s * m
```

until the oracle answers that the message s * m is PKCS compliant, which means there exists x for which s * m is in the interval [2 * B + x * N, 3 * B - 1 + x * N].
This means m is in [(2 * B + x * N) / s, (3 * B - 1 + x * N) / s]
The problem is we don't know x. However, we know candidates for x: these are all x for which interval [2 * B + x * N, 3 * B - 1 + x * N] has a non empty intersection with [s * 2 * B, s * (3 * B - 1)].
Thus, we know what are candidate intervals for m.

We need to find another s for which s * m is PKCS compliant. When we find it, we have a new set of candidate intervals for m. We compare intervals from both steps and if we find only one from the first step and one from the second step which have a non empty intersection, we have a new interval (a better estimation than [2 * B, 3 * B - 1]) where m lies. If there are more than one non empty intersections, we continue to search another s.

We continue finding s (in each step we double the s: s = 2 * s) and whenever a new compliant message is found, we narrow down the interval until it is of length 1 (this is message m). Note that after we narrowed down the interval for the first time (using the first two s), we are able to find new s values for which s * m is compliant much faster.

## Trimming

Improvements of the Bleichenbacher attack have been devised through years, but the most impressive is most likely the trimming technique invented by Bardou et all. [5].

The idea is that we find u and t for which u * m = x * t and GCD(u, t) = 1, where x is in some known range (m is plaintext that we are trying to get).

We know that m is in [2 * B, 3 * B] because it is PKCS compliant.
If we find such u and t and we for example know that x is also in the interval [2 * B, 3 * B], then we know that:

```  
t/u * 2 * B < m < t/u * 3 * B
```  

How do we find such u, t, x?
We simply try different u, t. For each (u, t) pair we calculate t_inv (inverse modulo n) and send u * t_inv * m to the oracle.
If the message is PKCS compliant, it means u * t_inv * m (mod n) is in the interval [2 * B, 3 * B].
We define x = u * t_inv * m (mod n). We know x is in [2 * B, 3 * B] when x is PKCS compliant.

If we take small enough u and t to u * m and x * t not exceed n, then we have:

```
u * m = x * t (mod n)
u * m = x * t
```

and we can narrow down the interval where m resides.

Note that m and x are smaller than 3 * B and B is much smaller than n (n is approximately B * 255 * 255), so we have many u, t for which u * m and x * t are smaller than n.
Also, we want to try at most a few thousand different (u, t) pairs otherwise there will be too many oracle calls required for finding the trimmer (u, t) and in this case it might be simply better to use the original Bleichenbacher algorithm.

# CBC-MAC (Matasano challenge 49)

CBC-MAC constructs MAC using a block cipher. It works the same as CBC mode encryption with some block cipher algorithm, but when the message is encrypted, only the last block is taken as MAC. 

Simple, but there are some caveats. For example - while when encrypting it is crucial to always use a different initialization vector (IV), it is absolutely necessary to always use the same IV with CBC-MAC. Otherwise, an attacker can arbitrarily change the first block of the message and still provide a valid MAC - the attacker needs to use a new IV for which xor(new_iv, modified_msg_first_block) = xor(iv, msg_first_block), in this case both, an original message and a modified message will have the same tag (and different IV):

```
from_id = 2
to_id = 3
amount = 1000
msg = "from=#{%s}&to=#{%s}&amount=#{%s}" % (from_id, to_id, amount)
key = "1" * 16
iv = b'\0' * BLOCK_SIZE
mac = get_mac(key, iv, msg)

to_id = 4
modified_msg = "from=#{%s}&to=#{%s}&amount=#{%s}" % (from_id, to_id, amount)
new_iv = strxor(strxor(msg[:16], iv), modified_msg[:16])
mac1 = get_mac(key, new_iv, modified_msg)

print(mac == mac1) # this is true
```

Also, length extension attack can be executed - let us say we have a message m and its tag t. If m is built out of m_1, m_2, ..., m_n blocks, we can construct a new message with a valid tag:

`m_1, m_2, ..., m_n||padding, m_1 xor t, m2, ..., m_n` is a message which has the same signature t as the original message `m_1, m_2, ..., m_n`.

# Do not mix MACs and hash functions (Matasano challenge 50)

You should not use MAC as a hash function or vice versa. MAC functions need to be existentially unforgeable under chosen message attack, while hash functions require collision resistance (see [security_definitions.md](https://github.com/miha-stopar/crypto-notes/blob/master/security_definitions.md)).

However, MACs can be built out of hash functions (like HMAC) and hash functions can be built out of block ciphers like AES - see Davies–Meyer construction in [crypto_constructs.md](https://github.com/miha-stopar/crypto-notes/blob/master/crypto_constructs.md).

Let's see how terribly wrong is to use CBC-MAC as a hash function.
The following Javascript code:

```
alert('MZA who was that?');
```

"hashes" to 296b8d7cb78a243dda4d0a61d33bbdd1 under CBC-MAC with a key "YELLOW SUBMARINE" and 0 as IV.
We can modify this Javascript code above so that it will alert some other value and still have the same "hash" as the original code.

If we have 5 blocks x1, x2, x3, x4, x5, it holds:

```
cbc-mac(x1 || x2 || x3 || xor(cbc-mac(x1 || x2 || x3), x4) || x5) = cbc-mac(x4 || x5)
```

The first code is 30 bytes long, we add "//" (which marks a beginning of a comment in Javascript) which makes it of two blocks length - these two are ours x1 and x2.

Then we take the following code

```
alert('Ayo, the Wu is back!');//
```

which are our x4 and x5. We want to modify this code to have the same "hash" as the original code.
For x3 we choose whatever we want, but it needs to be printable (valid Javascript code). For this reason we choose 7 random printable bytes, because this will be padded (PKCS7) into a printable bytes (the padding will be 9 bytes long and it will consist of 9xTAB, because TAB is 9 - we choose TAB, because it is printable).
However, xor(mac(x1 || x2 || x3), x4) might not be printable (in most cases it isn't), so we try different 7 random printable bytes, until we find a printable xor(mac(x1 || x2 || x3), x4).

The code we want is:

```
alert('Ayo, the Wu is back!');//   $)or\t\t\t\t\t\t\t\t\t8g4FJ-{Y*[gRH{tfas that?');
```

# Compression ratio side-channel attack (Matasano challenge 51)

Let's say a request is of the following format:

```
"""POST / HTTP/1.1
Host: hapless.com
Cookie: sessionid=TmV2ZXIgcmV2ZWFsIHRoZSBXdS1UYW5nIFNlY3JldCE=
Content-Length: %s
%s""" % (len(body), body)
```

If the request is compressed before the encryption, we can find out a cookie if we are able to trigger requests with different body content and observe into what what length the requests are compressed and encrypted.
We need to smartly choose body content and trigger the request. The header (with a secret token) is attached by the victim's browser.
First, we trigger the requests with the following body values:

```
Cookie: sessionid=A
Cookie: sessionid=B
Cookie: sessionid=C
Cookie: sessionid=D
...
Cookie: sessionid=T
```

By `Cookie: sessionid=T` we should see that the whole request has smaller length after it was compressed and encrypted. This is because the compression was better because the string `Cookie: sessionid=T` appeared two times in the request (otherwise, only a string `Cookie: sessionid=` is exploited for the compression).
However, we do not always see the difference in length because of encryption padding. But there are ways to resolve this - for example if the length of the request is the same for `Cookie: sessionid=T` and `Cookie: sessionid=A`, we check how many additional characters need to appended to both versions to get an additonal encryption block - we should see that more characters need to be appended to `Cookie: sessionid=T` to get an additional length because the compression is better here.

Once we know that T is the first character of a token, we repeat the process using the following body contents:

```
Cookie: sessionid=TA
Cookie: sessionid=TB
Cookie: sessionid=TC
Cookie: sessionid=TD
...
```

# Iterated multicollisions (Matasano challenge 52)

To construct a larger hash value sometimes two (or more) hashes are concatenated. For example by having two hash functions f and g, a new hash h is constructed as chash(m) = f(m) || g(m).

However, chash has no better collision resistance than f or g. Let us say the hash length (in bits) of f is n_1 and the hash length of g is n_2 (this means there are 2\*\*(n_1) possible outputs for f and 2\*\*(n_2) possible outputs for g).

Due to the birthday paradox, if we choose 1.2 * 2**(0.5 * n_1) messages, the probability that there will be at least one f collision in these messages, is greater than 0.5.

If for some reason one function (f or g) is weak and we can find a collision for it, we can quite easily find a collision for chash as well.

First we need to know that if k collisions of a single block length can be found ({m_1i, m_2i} for i=1,...,k such that hash(m_1i) = hash(m_2i)), 2\*\*k collisions can be easily constructed out of the original k collisions.

Merkle-Damgard hash functions are built from one-way compression functions (see [crypto_constructs.md](https://github.com/miha-stopar/crypto-notes/blob/master/crypto_constructs.md)), note that there is a fixed initial state. If h is compression function and iv is initial state, a one-block collision is:

f(m_10) = h(iv, m_10) = h(iv, m_11) = f(m_11) where m_0 and m_1 are two different blocks. 

Let us denote h1 = f(m_10) = f(m_11). As we assumed f is weak we can find such two blocks. We then search for two different blocks m_20, m_21 such that: h(h1, m_20) = h(h1, m_21).

Let us denote this value as h2 = h(h1, m_20) = h(h1, m_21).

Now it holds: f(m_10 || m_20) = f(m_10 || m_21) = f(m_11 || m_20) = f(m_11 || m_21), thus we have four messages of block length two which all produce the same hash. This is due to:

```
f(m_10 || m_20) = h(h(iv, m_10), m_20) = h(h1, m_20) = h2
f(m_11 || m_20) = h(h(iv, m_11), m_20) = h(h1, m_20) = h2
f(m_10 || m_21) = h(h(iv, m_10), m_21) = h(h1, m_21) = h2
f(m_11 || m_21) = h(h(iv, m_11), m_21) = h(h1, m_21) = h2
```

Thus, if we want multicollision of block-length k, at each step we produce a collision of block-length 1 where we use temporary hash as initial state.

We then simply concatenate k blocks together, where for each block i we can choose either m_1i or m_2i (that is 2\*\*k different k-block messages). All these 2\*\*k messages have the same hash due to the Merkle-Damgard construction.

Now that we have a large pool of collisions for a weaker hash function, there is a good chance that among there is a collision for a stronger hash function as well (this is why we need a lot of collisions for a weaker hash).
If there is a collision for a stronger hash, let's say g(m1) = g(m2), then f(m1) || g(m1) = f(m2) || g(m2) and we have a collision for chash.

# Expandable messages attack (Matasano challenge 53)

The attack is described in [4].

Cryptographic hash function needs to be resistant to second preimage attacks. That means if you have x1 where hash(x1) = y, you should not be able to find x2 such that hash(x2) = y.

Now let us say we have a very long message (n = 2**k + k - 1 blocks): 

```
msg = block_0 || block_1 || ... || block_(n-1).
```

We calculate all intermediate hashes: iv=H_0, H_1, ..., H_n (see Merkle-Damgard in [crypto_constructs.md](https://github.com/miha-stopar/crypto-notes/blob/master/crypto_constructs.md)). Let us say c is a compression function used in Merkle-Damgard construction.

H_k is an intermediate hash after processing k blocks.

Let's assume we can find a collision of some block with some of the intermediate hashes after k blocks (if k is large, brute-force attack becomes feasible):

```
c(H_k, connecting_block) = H_l where l > k.
```

Note that we used H_k as the initial state. This means that if no length padding is applied:

```
hash(block_0 || ... || block_(k-1) || connecting_block) = hash(block_0 || .. || block_(l-1))
```

and also:

```
hash(block_0 || ... || block_(k-1) || connecting_block || block_l || ... || block_(n-1)) = hash(msg)
```

However, if length padding is used this would not hold. But we can expand first k blocks into l-1 blocks which would result in the same hash as H_k. This is where expandable messages come into play. We need l-1 blocks such that:

```
H_k = hash(new_block_0 || ... || new_block_(l-2))
```

Let us see how we can replace block_i with 2**i + 1 blocks as follows:

We choose some dummy_message of block lenth 2**i and calculate iv_d = hash(dummy_message) with initial state of H_i (and without appending padding to the dummy_message). Then we search for a collision such that (let the second parameter in hash function be the initial state):

```
hash(some_block, iv_d) = H_(i+1) = hash(block_0 || ... || block_(i-1) || dummy_message || some_block)
```

This also means:

```
hash(block_0 || ... || block_(i-1) || dummy_message || some_block || block_(i+1) || ... || block_(k-1)) = H_k
```

And also (if no length padding):

```
hash(block_0 || ... || block_(i-1) || dummy_message || some_block || block_(i+1) || ... || block_(k-1) || connecting_block || block_l || ... || block_(n-1)) = hash(msg)
```

We expanded first k blocks by 2\*\*i blocks. If we do the same for each i < k, we can expand first k blocks by any number between 1 and 2\*\*k - 1 (by properly combining which blocks to expand and which not) and thus obtain a message of the same length as msg.

# Herding attack (Matasano challenge 54)

Let's say we want to bet on a result of some basketball game. For some reason we are concerned that publishing our prediction will affect the outcome of a game. So we publish only a hash of our prediction message.

It turns out [3] that if the hash function is used for which many collisions can be found, we can construct after the game ends a message which will contain the proper result (which we now know) and will have the same hash as we submitted before the game started.

For this we need to prepare a diamond structure [3] first.
We generate 2^k messages with different hashes. Parameter k needs to be large, but for demonstration purposes let's say k = 3.
So we have m_01, m_02,..., m_08 with hashes h_01, h_02,...,h_08. We now choose m_11 and m_12 such that

```
hash(m_11, h_01) = hash(m_12, h_02) where the second argument is initial value
```

this means if we concatenate m_01 with m_11 and m_02 with m_12 we get the same hash (default initial value): 

```
hash(m_01 || m_11) = hash(m_02 || m_12)
```

We choose m_13, m_14,...,m_18 analogously (hash(m_03 || m_13) = hash(m_04 || m_14) ...).

```
m_01   | m_11 (h_11)
m_02   | m_12 (h_11)
m_03   | m_13 (h_12) 
m_04   | m_14 (h_12)
m_05   | m_15 (h_13)
m_06   | m_16 (h_13)
m_07   | m_17 (h_14)
m_08   | m_18 (h_14)
```

The value in brackets is hash of the concatenation of the messages from this row.
We then choose m_21 and m_22 such that:

```
h_21 := hash(m_21, h_11) = hash(m_22, h_12) where the second argument is initial value
```

This means:

```
hash(m_01 || m_11 || m_21) = hash(m_02 || m_12 || m_21) = hash(m_03 || m_13 || m_22) = hash(m_04 || m_14 || m_22)
```

Analogously for m_23 and m_24. We now have:

```
m_01   | m_11 (h_11)  | 
m_02   | m_12 (h_11)  | m_21 (h_21)
m_03   | m_13 (h_12)  |
m_04   | m_14 (h_12)  | m_22 (h_21)
m_05   | m_15 (h_13)  |
m_06   | m_16 (h_13)  | m_23 (h_22)
m_07   | m_17 (h_14)  |
m_08   | m_18 (h_14)  | m_24 (h_22)
```

We continue so that:

```
m_01   | m_11 (h_11)  | 
m_02   | m_12 (h_11)  | m_21 (h_21)
m_03   | m_13 (h_12)  |
m_04   | m_14 (h_12)  | m_22 (h_21)  | m_31 (h_31)
m_05   | m_15 (h_13)  |
m_06   | m_16 (h_13)  | m_23 (h_22)
m_07   | m_17 (h_14)  |
m_08   | m_18 (h_14)  | m_24 (h_22)  | m_32 (h_31)
```

We submit h_31 as our prediction. Once the result of the game is known, let's say 17:11, we search for a link_block such that:

```
m = "it will be 17:11"
hash(m || link_block) = one of the intermediate hashes in our diamond structure
```

Let's say we found a link_block such that:

```
hash(m || link_block) = h_14
```

Now:

```
hash(m || link_block || m_24 || m_32) = h_31
```

Thus, we found a message which contains the proper prediction and has the hash that we submitted (note that there might be a problem with message length and padding - in this case expandable messages can be used, see the previous challenge).

# MD4 collisions (Matasano challenge 55)

MD4 compresses any arbitrary bit-length message into a 128-bit hash value. It first pads the message into a message with a length that is a multiple of 512 bits. For each 512-bit message block, it compresses it into a 128-bit hash value using a compression function (the output of one block is used in processing of the next one).

Let us see how compression function works on a 512 bit-length block m. Block m is split into 16 words (each word is of 32-bit length) - m_0, m_1, ..., m_15.

Compression function has three rounds. Each round uses a different nonlinear Boolean function: F, G, and H. The operations of these three functions are all bitwise. There is a state h which gets updated in each step of each round (each round has 16 steps - in each step one of the 16 words is used). The initial value for h is fixed. State h consists of 4 words (in each step one of the words is updated).

First round iterates over block words as follows:

```
h[0] = leftrotate((h[0] + F(h[1], h[2], h[3]) + m_0) % 2**32, 3) 
h[3] = leftrotate((h[3] + F(h[0], h[1], h[2]) + m_1) % 2**32, 7) 
h[2] = leftrotate((h[2] + F(h[3], h[0], h[1]) + m_2) % 2**32, 11) 
h[1] = leftrotate((h[1] + F(h[2], h[3], h[0]) + m_3) % 2**32, 19) 
h[0] = leftrotate((h[0] + F(h[1], h[2], h[3]) + m_4) % 2**32, 3) 
...
```
It turns out (see Wang [2]) that if message msg is slightly changed in some particular way to msg1, both messages have the same hash. In this case there is a minor difference in h state which is propagated through all the steps and rounds and eventually becomes 0 at step 21 - at steps 36 and 37 some difference appear again but is zeroed again (however, not in all cases).

Message msg1 has to be constructed from msg as follows (we consider only messages of 512 bit length now - one block):

```
msg1_1 = msg1_1 + 2^31
msg1_2 = msg1_2 + 2^31 - 2^28
msg1_12 = msg1_12 - 2^16
```

where msg1_1 is the first word of a block, msg1_2 the second...

Let us denote the temporary state of msg1 (of the modified message) as hh. 

We can see that after first step h[0] is the same as hh[0] (because at this point h[1]=hh[1], h[2]=hh[2], h[3]=hh[3], and msg_0=msg1_0).

After second step h[3] and hh[3] are different because msg_1 != msg1_1. The difference is 2^6 - this is 2^31 rotated to the left by 7.

In the third step we have for the first time a case where F arguments differs for msg and msg1 (actually only h[3] differs). Wang attack says that a difference after the third step (difference between h[2] and hh[2]) should be 2^10 - 2^7. We see that this is exactly the difference in msg_2 and msg1_2 leftrotated by 11. Thus we need: F(h[3], h[0], h[1]) = F(hh[3], hh[0], hh[1]).

Note that h[0]=hh[0], h[1]=hh[1], and hh[3]-h[3] = 2^6. That means h[3] and hh[3] differs only in a bit at position 7 (if the least significant bit is at position 0). As F is bitwise (as well G and H) we have to check what needs to be done to make F(h[3], h[0], h[1]) = F(hh[3], hh[0], hh[1]) in a bit at position 7 (for other bits it is already the same).

F is defined as:

```
F(x, y, z) = (x and y) or (not x and z)
```

In our case we have:

```
F(h[3], h[0], h[1]) = (h[3] and h[0]) or (not h[3] and h[1])
F(hh[3], hh[0], hh[1]) = (hh[3] and hh[0]) or (not hh[3] and hh[1])
```

We know:

```
F(hh[3], hh[0], hh[1]) = F(hh[3], h[0], h[1]) = (hh[3] and h[0]) or (not hh[3] and h[1])
```

This is the same as F(h[3], h[0], h[1]) when h[0]=h[1] (at position 7). This is the first of the conditions for msg for this attack to work (if the conditions are not met, the modifications to get msg1 won't help; note that messages which for which the conditions hold are rare).

In the fourth step we calculate h[1]:

```
h[1] = leftrotate((h[1] + F(h[2], h[3], h[0]) + m_3) % 2**32, 19) 
```

We need hh[1]-h[1] = 2^25 (the required differences at each step are called differential path). We have: h[0]=hh[0], h[1]=hh[1], h[2]!=hh[2], and h[3]!=hh[3].

Thus, we need F(hh[2], hh[3], hh[0]) - F(h[2], h[3], h[0]) = 2^6 (if we leftrotate 2^6 by 19 we get 2^25).

This means F(hh[2], hh[3], hh[0]) and F(h[2], h[3], h[0]) need to the same in all bits except in bit at position 7. We know that h[0]=hh[0] in all bits; h[2] and hh[2] are different in position 8 and 11; h[3] and hh[3] are different in position 7.

We need to check positions 7, 8, 11.

At position 7 we have difference in h[3] and hh[3]. If we want to retain the difference h[2] needs to be 1 at this position because F(x, y, z) = F(x, y1, z) when x=1 (here y!=y1).

At position 8 we have difference in h[2] and hh[2]. We want to get rid of this difference. If we analyze F we see that y needs to be the same as z to achieve this and thus h[3] needs to be the same as h[0] at position 8.

At position 11 we have difference in h[2] and hh[3]. Analogously as above h[3] and h[0] needs to be the same at position 11.

This way the required differences (differential path) can be provided for the first round (first 16 steps) - if the message msg meets the conditions. However, the differential path needs to be followed in the second and third round as well to eventually get a collision. 
And this is harder, because once we are in a second round we cannot freely change the message blocks to get the required differences - if we change a block while in a second round, we ruin what we have done in the first round (we might break some differences in the first round).

However, there are some tricks that we can still apply (multi-step modifications). We can change the block as we need to get a required difference (in a step in round two) and then try to minimize the effect in round one. Obviously, if we change for example a block 1 to get a required difference in round 2, the state hh after the second step in the first round (where block 1 is used) will be changed (hh[3] is changed).
However, the required difference hh[3]-h[3] may still be there, but the problem is that hh[3] is used also used in the next four steps (to compute hh[2], hh[1], hh[0], and hh[3]).
But now we can change block 2 to neutralize the change when computing hh[2] in the step three (that means hh[2] is the same as it was prior the modification of block 1). The same we do with blocks 3, 4, and 5 - after the sixth step all four elements of the hh are exactly the same as without the modified block 1.
The only difference stayed in hh[3] after the second step, but this does not always ruin the differential path. However, sometimes it does and in these cases this method doesn't work.

There are other tricks - for example changing the late blocks, like the last block (15) where neutralizing is not needed, since we didn't break the differential path in round 1. However, in most cases this becomes even more complex than with the method above.
Bottom line, not for all messages collisions can be found. Wang [1] reports that 4 to 64 randomly selected messages are needed to find a colision, however my implementation needs a few thousand randomly selected message to find one - probably I didn't use all the multi-step modification tricks I could.

# RC4 biases (Matasano challenge 56)

RC4 is a stream cipher - it takes a key and generates a stream of bytes which is then used to xor the plaintext and thus to encrypt it.
The stream of bytes needs to be random, but in RC4 it is not. There are small biases (see [1]), for example the byte at position 16 tends to be 240, the byte at position 32 tends to be 224 ... 
There are also multi-byte biases, but already with the two single-byte biases mentioned above are pretty powerful. 

For example for decryption of a secret cookie which is included into some web requests. Let's assume the cookie is at the beginning of the request (which is not in reality), now if we can somehow manage to execute many of these requests (let's say 2^22) and monitor the traffic, we can observe to which byte the 16-th byte of a token is encrypted. 
Due to the bias, this byte is supposed to be xored with byte 240 in slightly more cases than with other bytes. So if we take the byte at 16-th position that appeared in the most cases, and we xor it with 240, we should get an actual 16-th byte of a secret token.
For 15-th byte we need to somehow prepend a byte to a request and then repeat the process above.



[1] AlFardan, N., Bernstein, D., Paterson, K., Poettering, B., Schuldt,
J.: On the Security of RC4 in TLS. http://www.isg.rhul.ac.uk/tls/
index.html (2013). Accessed Apr. 2013

[2] X. Wang, X. Lai, D. Feng, H. Chen, and X. Yu. Cryptanalysis of the hash functions MD4
and RIPEMD. In R. Cramer, editor, Advances in Cryptology – EUROCRYPT 2005,
volume 3494 of Lecture Notes in Computer Science. Springer-Verlag, Berlin, Germany,
2005.

[3] T. Kohno and J. Kelsey, Herding hash functions and the Nostradamus attack,
Advances in Cryptology – Eurocrypt 2006 (S. Vaudenay, ed.), LNCS, no. 4004,
Springer-Verlag, 2006, pp. 222–232.

[4] J. Kelsey and B. Schneier. Second preimages on n-bit hash functions for much less
than 2n work. In R. Cramer, editor, Advances in Cryptology – EUROCRYPT 2005,
volume 3494 of Lecture Notes in Computer Science, pages 474–490. Springer-Verlag,
Berlin, Germany, 2005.

[5] R. Bardou, R. Focardi, Y. Kawamoto, L. Simionato, G. Steel, and J.-K. Tsay. Efficient padding oracle attacks on cryptographic hardware. In CRYPTO, pages 608–625, 2012
