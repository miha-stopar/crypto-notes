Convert a Unicode number into a character (start node.js console):
> String.fromCharCode(22)
'\u0016'

> String.fromCharCode(65)
'A'

Convert character to ASCII:
> 'h'.charCodeAt(0)

Convert ASCII value to character:
> String.fromCharCode(104)
'h'

Javascript interprets a number as hexadecimal if there is 0x prepended:
> 0xFF
255

# Floating point
	
Javascript does not natively support arbitrary-precision (bignum) arithmetic which is needed for cryptography (for example to prevent overflows when doing arithmetic with integers having hundreds of digits is needed), but there are many Javascript BigInteger libraries available.

Javascript natively is not calculating arbitrary large integeres precisely. It is accurate up to 15 digits. See the example below:

> x = 999999999999999
999999999999999
> x = 9999999999999999
10000000000000000

Also:
> Math.pow(2,53)
9007199254740992
> Math.pow(2,53)+1
9007199254740992
> Math.pow(2,53)+2
9007199254740994
> Math.pow(2,53)+3
9007199254740996
> Math.pow(2,53)+4
9007199254740996

This is because the numbers are stored in 64 bits: first 52 bits are used to store the fraction (mantissa), the exponent is stored in bits from 52 to 62, the last bit is used to store the sign. Thus the operations on numbers bigger than Math.pow(2, 53) do not execute necessarily correctly.

Also, when operating with decimals, Javascript is not precise. Open node.js console and type:
> 0.1+0.2
0.30000000000000004

You get something that is slightly different than expected. This is because floating point numbers are represented in hardware as base 2. And most decimal fractions cannot be represented exactly as binary fractions.

You can represent for example 0.125, which is 0.001 (2 to the power of -3). However most of the decimals are not of this format (or sum of number of this format).

Just as a curiosity: while many programming languages define different types of numbers (integeres, short, long...), Javascript numbers are always stored as (double precision) floating point numbers.

# Bitwise operators

Bitwise operators converts the number to a 32-bit signed integer and thus discards any fractions and higher-placed bits than 32.

(Math.pow(2,32) + 10)>>32
10







