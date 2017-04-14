#AES 

AES operates on 4x4 matrix of bytes which is named a state (state can contain 16 bytes which is the same as 128 bits - which is AES block size). Likewise, the cipher key is contained in a similar matrix with four rows and 4, 6, or 8 columns (depends on the key size). Each element in the matrix can be considered as an element from GF(2^8). More about GFs can be found in [number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md).

AES consists of a number of transformation rounds that convert the plaintext into the ciphertext. The number of transformation rounds depend on the key size (128: 10 rounds, 192: 12 rounds, 256: 14 rounds). All rounds except the last one are identical.

Each round consists of a single-byte based substition (SubBytes), a row-wise permutation (ShiftRows), a column-wise mixing (MixColumns), and the addition of the round key (AddRoundKey).

First, the input block (16 bytes) is put into a 4x4 matrix so that the first four bytes are in the first column, the next four bytes in the second column... Similarly the key is arranged in the matrix.

While the matrix is named a state, there is also a notion of word, which is a sequence of four bytes (a column is a word, as well is a row).

In what follows the AES 128 (128-bit key) is briefly explained. The four columns (four words) of it are expanded into 44 words key schedule - 4 words for each round.

## Cipher

At the start of the cipher, the input is arranged into the State as described above. First, an initial Round Key addition takes place, then the State is transformed by applying the round function 10, 12, or 14 times. The final round is slightly different than the rest - it does not include the MixColumns step.

The round is parameterized with a key schedule.

## SubBytes

SubBytes uses 16x16 substitution table to find a replacement byte for a given byte in the input state. Substitution table is created on the basis of GF(2^8) and bit scrambling. It is up to this step to reduce the correlation between the input bits and the output bits. The bit scrambling is needed to ensure that the substitution cannot be evaluated as a simple function.

Substition table (S-box) is constructed by composing two transformations. First we prepare the initial table where value in each cell is the joined value of row index and column index:

![sbox](https://raw.github.com/miha-stopar/crypto-notes/master/img/sbox1.png)

Then we replace the value in each cell by its multiplicative inverse in GF(2^8), using the irreducible polynomial x^8 + x^4 + x^3 + x + 1. The value 0x00 is replaced by itself (it does not have a multiplicative inverse). After this step, bit scrambling is applied on the value on each cell. The value in a cell should now be thought as a vector, like 10010001 (coordinates b_7, b_6,...,b_0). Bit scrambling is done using the affine transformation:

![sbox](https://raw.github.com/miha-stopar/crypto-notes/master/img/sbox3.png)

Once the bit scrambling is applied, the result is S-box:

![sbox](https://raw.github.com/miha-stopar/crypto-notes/master/img/sbox2.png)

When applying SubBytes for a particular value, say 53 (which in decimal is 5 * 16 + 3 = 83), simply the value in the matrix at row index 5 and column index 3 is taken, which is ed (237 in decimal).

### Attacks on S-boxes

The only nonlinear part of the most block ciphers are S-boxes. Thus, it is crucial that S-Box cannot be approximated by a linear function - there should not exist a linear relationship between a subset of the plaintext bits and a subset of the ciphertext bits (the relationship should hold with probability 0.5 + epsilon, where epsilon is non-negligible value, called bias), as this would make the cipher vulnerable to known-plaintext attacks (for more about type of attacks check [security_definitions.md](https://github.com/miha-stopar/crypto-notes/blob/master/security_definitions.md)). 

For example there was a linearity found in DES with a bias greater than 2^(-21). It was shown that with 2^(43) given plaintext/ciphertext pairs, the attack succeeds with probability of 85%. However, having 2^(43) known pairs is not realistic as this correspond to the 64 terabytes of plaintext data, but still it shows how important the nonlinearity of S-boxes (and of course consequently of the whole cipher) is.

## ShiftRows

ShiftFows does not change the first row, shifts the second row by one byte to the left, the third row by two bytes, and shifts the last row by three bytes (it scrambles the byte order inside the 128-bit block):

![sbox](https://raw.github.com/miha-stopar/crypto-notes/master/img/sbox4.png)

## MixColumns

MixColumns further scrambles the byte order. Both, ShiftRows and MixColumns take care that each bit of the ciphertext depends on every bit of the plaintext after all rounds are finished.

MixColumns takes the column and each byte in it replaces by two times that byte, plus three times the next byte in the column, plus the byte that follows, plus the byte that follows (it goes in cycles - the byte that follows the fourth element is again the first one). For example the new element at position 00 would be:
0x02 * s_00 + 0x03 * s_10 + s_20 + s_30


This is obviously multiplication and addition in GF(2^8) (addition is xor). That means that multiplying by two is multiplying by the polynomial that corresponds to 00000010, and multiplying by three is multiplying by the polynomial that corresponds to 00000011.

If for example we have a column [56, 169, 150, 68], we calculate:

```
0x02 * 56 + 0x03 * 169 + 150 + 68 = 112 + 224 + 150 + 68 = 66
```

So the first element in a column is replaced by 66. Note that:

```
0x03 * 169 = 00000011 * 10101001 = (x + 1) * (x^7 + x^5 + x^3 + 1) = (x^8 + x^6 + x^4 + x) + (x^7 + x^5 + x^3 + 1) = 
= x^6 + x^4 + x + (x^4 + x^3 + x + 1) + x^7 + x^5 + x^3 + 1 = x^7 + x^6 + x^5 + 2 * x^4 + 2 * x^3 + 2 * x + 2 = x^7 + x^6 + x^5 = 192 + 64 + 32 = 224
```

For information how to multiply in GF(256) see [number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md).

## AddRoundKey

AddRoundKey takes the round key (derived from the original 128-bit key) and xors it with the state.

### Key expansion

First, the 16 bytes of the original key are arranged in the matrix in the same way it is done for the input block. 

![key expansion](https://raw.github.com/miha-stopar/crypto-notes/master/img/key_expansion1.png)

The first column (the first four bytes of the key) constitute the word w_0, the second w_1...:

The key expansion algorithm expands the four words to 44 words. The four words [w_0, w_1, w_2, w_3] are bitwise xored with the input block before the rounds start. The remaining 40 words are used in the 10 rounds (for AES-128).

The round key is derived from the key in previous round. If we have [w_i, w_(i+1), w_(i+2), w_(i+3)], we can generate [w_(i+4), w_(i+5), w_(i+6), w_(i+7)]. This key will be used in (i/4)-th round.

Calculating w_(i+5), w_(i+6), w_(i+7) goes simply like:

```
w_(i+5) = w_(i+4) + w_(i+1)
w_(i+6) = w_(i+5) + w_(i+2)
w_(i+7) = w_(i+6) + w_(i+3)
```

Calculating w_(i+4) is done as:

```
w_(i+4) = w_i + g(w_(i+3))

```

Function g consists of three steps:

 * one-byte circular rotation of 4-byte word 
 * byte substitution of the 4 bytes of the word returned by previos step (circular rotation) using S-box
 * xor the bytes from the previous step with a round constant (see below)

#### Round constant

Round constant is a word with the three rightmost bytes being zero. The round constant for the i-th round is denoted Rcon[i]:

```
Rcon[i] = [RC[i], 0x00, 0x00, 0x00]
```

RC is defined as:

 
```
RC[1] = 0x01
RC[j] = 0x02 x RC[j-1]
```

The purpose of the round constants is to break any symmetries that may have appeared in the other steps of the key expansion algorithm.





