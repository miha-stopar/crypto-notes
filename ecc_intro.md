# Intorduction to the elliptic curve cryptography

## Weierstrass equation

Elliptic curve E over finite field K is defined by a Weierstrass equation:


```
y^2 + a_1*x*y + a_3*y = x^3 + a_2*x^2 + a_4*x + a_6 
```

where a_1, a_2, a_3, a_4, a_6 are from K and the discrimant of E is not 0.

Before giving a definition of a discriminant, let us see how the equation changes if we change a variable `y` into `2*y + a_1*x + a_3` and assuming that field characteristic is not 2:

```
y^2 = 4*x^3 + b_2*x^2 + 2*b_4*x + b_6 
```

where:

```
b_2 = a_1^2 + 4*a_2
b_4 = 2*a_4 + a_1*a_3
b_6 = a_3^2 + 4*a_6
```

We also define:

```
b_8 = a_1^2*a_6 + 4*a_2*a_6 -a_1*a_3*a_4 + a_2*a_3^2 - a_4^2
```

A discriminant is defined as:

```
discr = -b_2^2*b_8 - 8*b_4^3 - 27*b_4*3 + 9*b_2*b_4*b_6
```

If the field characteristic is neither 2 or 3, the equation can be transformed into (simplified Weierstrass equation):

```
y^2 = x^3 + a*x + b

```

Then, a discriminant is `-16(4*a^3 + 27*b^2)`.

A nonzero discriminant algebraically means that a right-hand side of the equation has three distinct roots. Geometrically it means that the elliptic curve is smooth - there are no points at which the curve has two or more distinct tangent lines.

## Elliptic curve cryptography

For cryptographic purposes an elliptic curve consists of the points satisfying the above equation (simplified Weierstrass) along with a point at infinity.

Let us see an elliptic curve in sage:

```
sage: F = EllipticCurve(GF(19), [-7, 10])
sage: F
Elliptic Curve defined by y^2 = x^3 + 12*x + 10 over Finite Field of size 19
```
Note that -7 is the same as 12 in GF(19).

```
sage: F.count_points()
24
sage: F.points()
[(0 : 1 : 0), (1 : 2 : 1), (1 : 17 : 1), (2 : 2 : 1), (2 : 17 : 1), (3 : 4 : 1), (3 : 15 : 1), (5 : 9 : 1), (5 : 10 : 1), (7 : 0 : 1), (9 : 7 : 1), (9 : 12 : 1), (10 : 3 : 1), (10 : 16 : 1), (12 : 1 : 1), (12 : 18 : 1), (13 : 8 : 1), (13 : 11 : 1), (16 : 2 : 1), (16 : 17 : 1), (17 : 4 : 1), (17 : 15 : 1), (18 : 4 : 1), (18 : 15 : 1)]
```

The third coordinate specifies whether the point is at infinity (0 means infinity). We need to have a point at infinity as identity element.

You can plot the curve in sage as well:

```
sage: P = plot(F, rgbcolor=(0,0,1))
sage: P
Launched png viewer for Graphics object consisting of 1 graphics primitive
```

![ECC](https://raw.github.com/miha-stopar/crypto-notes/master/img/ecc1.png)

An example elliptic curve is given (all points drawn, not only from K) here (taken from [Wikipedia](https://en.wikipedia.org/wiki/Elliptic_curve)):

![ECC addition](https://raw.github.com/miha-stopar/crypto-notes/master/img-wikipedia/ecc-addition.png)

The picture above also demonstrates how addition of two points is defined for elliptic curves - a line needs to be drawn from P to Q, the intersection of this line with the curve is then relected over the x axis to get -R = P + Q (picture shows only R). This is called chord-and-tangent rule for adding points. Such addition turns elliptic curve over a finite field into an abelian group where infinity (denoted as 0) serves as an identity element.

Algebraically, the addition can be calculated as follows. Let us say that P = (x1, y1), Q = (x2, y2) and R = (x3, y3). We need to find x3 and y3. R will then be (x3, -y3).

We know that R lies on the curve, so for (x1, y1), (x2, y2), and (x3, y3) it holds:

```
y^2 = x^3 + a*x + b mod p
```

And R is on the line through P and Q, let us say this line is:

```
y = k*x + n
```

We can easily calcalute k and n (we actually will not need n):

```
k = (y2-y1)/(x2-x1) 
```

By solving these two equations we get (we omit modulo p): 

```
(k*x + n)^2  = x^3 + a*x + b
k^2*x^2 + 2*k*x*n + n^2 = x^3 + a*x + b
x^3 - k^2*x^2 + a*x - 2*k*x*n - n^2 + b = 0
x^3 - k^2*x^2 + (a - 2*k*n)*x - n^2 + b = 0
```

This cubic polynomial has three roots (x1, x2, x3) and we know two of them. Thus

```
x^3 - k^2*x^2 + (a - 2*k*n)*x - n^2 + b = (x-x1)*(x-x2)*(x-x3)
x^3 - k^2*x^2 + (a - 2*k*n)*x - n^2 + b = (x-x1)*(x-x2)*(x-x3)
x^3 - k^2*x^2 + (a - 2*k*n)*x - n^2 + b = x^3 - (x1 + x2 + x3)*x^2 + (x1*x2 + x1*x3 + x2*x3)*x - x1*x2*x3
```

If we compare coefficients by x^2 we get:

```
k^2 = x1 + x2 + x3
x3 = k^2 - x1 - x2
```

We can calculate y3 easily:

```
y3 = y1 + k*(x3-x1)
```

Thus, we obtained coordinates for R:

```
x3 = k^2 - x1 - x2
y3 = y1 + k*(x3-x1)
```

For -R = P + Q:

```
x3 = k^2 - x1 - x2
y3 = - y1 - k*(x3-x1)

```

However, there are some special cases (images from 2 to 4), for example where P = Q, which are covered in images 2 - 4.

For P = Q, the slope of tangent line at P is computed. Let us say E(x, y) = x^3 + a * x + b - y^2 and P = (x1, y1). The slope is then:

```
lambda = [∂E/∂X / ∂E/∂Y](x1, y1) = (3 * x1^2 + a) / (2 * y1)
```

The point is:

```
x3 = lambda^2 - 2 * x1
y3 = (x1 - x3) * lambda - y1
```

## How to speed up ECC operations

### Reduction

Reduction z mod p can be expensive (note that in the operations above for adding two points there is always a reduction mod p). To speed up the reduction NIST recommeded certain primes (to be used as prime p) which enable fast reduction.

NIST primes are the following:

 * p_192 = 2^192 - 2^64 - 1
 * p_224 = 2^224 - 2^96 + 1
 * p_256 = 2^256 - 2^224 + 2^192 + 2^96 - 1
 * p_384 = 2^384 - 2^128 -2^96 + 2^32 - 1
 * p_521 = 2^521 - 1

As can be seen these primes can all be written as a sum or difference of a small number of powers of 2. Also (except for p_521) the powers in these expressions are all multiples of 32. These properties enable efficient reduction and are especially convenient for machines with wordsize 32.

For example if we take p = p_192 = 2^192 - 2^64 - 1 and a positive integer c < p^2, we can write c in base 2^64 presentation as:

```
c = c_5 * (2^64)^5 + c_4 * (2^64)^4 + c_3 * (2^64)^3 + c_2 * (2^64)^2 + c_1 * 2^64 + c_0
```

where c_i < 2^64.

When in Z_p we know (as p = 2^192 - 2^64 - 1):

 * 2^192 = 2^64 + 1
 * 2^256 = 2^128 + 2^64 (multiplying the first equation by 2^64)
 * 2^320 = 2^128 + 2^64 + 1 (multiplying the second equation by 2^64 and applying the first equation to replace 2^192)

Thus, c (in Z_p) is:

```
c = c_5 * 2^128 + c_5 * 2^64 + c_5 + c_4 * 2^128 + c_4 * 2^64 + c_3 * 2^64 + c_3 + c_2 * 2^128 + c_1 * 2^64 + c_0
```

Now we can repeatedly subtract p until we get the result less than p.

#### Example: calculations in P224 (this is how golang's p224.go is implemented)

As we saw: `p_224 = 2^224 - 2^96 + 1`. We can write each element from Z_p (where p = p_224) as follows (because (2^28)^8 = 2^224):

```
a = a_7 * (2^28)^7 + a_6 * (2^28)^6 + ... + a_1 * 2^28 + a_0
```

Number used in cryptography are too large to fit into a single CPU register. This is why it is more efficient to use for example base 2^28. Coefficients a_i are called limbs. Base 2^32 could be used to, but then there are problems with overflowing. 

##### Addition

Let's have a look at addition. If `b = b_7 * (2^28)^7 + b_6 * (2^28)^6 + ... + b_1 * 2^28 + b_0`, then a + b is:

```
a + b = (a_7 + b_7) * (2^28)^7 + (a_6 + b_6) * (2^28)^6 + ... + (a_1 + b_1) * 2^28 + (a_0 + b_0)

```

While a_i and b_i are all smaller than 2^28, a_i+b_i are smaller than 2^29. This means limbs of a+b are not necessarily smaller than 2^28 and we don't have a representation in base 2^28 anymore. However, the limbs still fit into a processor register and we can efficiently work with this presentation. When base 2^28 is needed, there is some work to convert a+b back to have all limbs under 2^28.

##### Subtraction

We could calculate a-b as:

```
a - b = (a_7 - b_7) * (2^28)^7 + (a_6 - b_6) * (2^28)^6 + ... + (a_1 - b_1) * 2^28 + (a_0 - b_0)

```

However, a_i - b_i might be negative and we would like to simplify this and have only positive limbs.

Thus, we choose a multiplication of p224, so that all limbs are very close to 2^31. It turns out 8*p224 is what we are looking for:

```
mulZero = 8 * p224 = (2^31 + 2^3) * (2^28)^7 + (2^31 - 2^3) * (2^28)^6 + (2^31 - 2^3) * (2^28)^5 + (2^31 - 2^15 - 2^3) * (2^28)^4 + (2^31 - 2^3) * (2^28)^3 + (2^31 - 2^3) * (2^28)^2 + (2^31 - 2^3) * (2^28)^1 + (2^31 - 2^3)

```

As we are in field Z_p224, mulZero is a multiplication of zero and thus zero itself. We can do:

```
a - b = (a_7 - b_7 + mulZero_7) * (2^28)^7 + (a_6 - b_6 + mulZero_6) * (2^28)^6 + ... + (a_1 - b_1 + mulZero_1) * 2^28 + (a_0 - b_0 + mulZero_0)

```

Now all limbs are positive and smaller than 2^32.

##### Inversion

We could extended Euclidean algorithm (see [number_theory.md](https://github.com/miha-stopar/crypto-notes/blob/master/number_theory.md)) to invert an element a, however this is not constant time. Instead, we calculate it as:

```
a^(-1) = a^(p224 - 2)
```

This works because a^p = a (mod p) when p is primer (Fermat's little theorem).




### Transformation between coordinates

To calculate k = (y2-y1)/(x2-x1) we need to calculate modular inverse of (x2-x1) which is then multiplied by (y2-y1). Modular inversion is computationally expensive and for this reason the affine coordinates are usually transformed into projective or Jacobian coordinates.

Conversion to Jacobian coordinates can be done as follows. Let's first introduce an equivalence relation ~ on the set K^3\{(0,0,0)}:

```
(x1, y1, z1) ~ (x2, y2, z2) if x1 = lambda^c * x2, y1 = lambda^d * y2, z1 = lambda * z2 for some lambda from K*
```

The equivalence class is called a projective point.

We have P(x, y) = y^2 - x^3 - a * x^2 - b. We are interested in (x, y) such that P(x, y) = 0. 

Now we would like to have P(x, y, z) such that if P(x1, y1, z1) = 0, then P(x2, y2, z2) = 0 for all (x2, y2, z2) that are equivalent to (x1, y1, z1). And also that P(x, y, 1) = 0 if P(x, y) = 0.

We get this by (x, y) -> (x/(z^c), y/(z^d)). For c = 1, d = 1:

```
y^2 * z = x^3 + a * x^2 * z + b * z^3
```

Example: y^2 = x^3 + a * x + b and c = 2, d = 3.
After transformation we get:

```
y^2 = x^3 + a * x * z^4 + b * z^6
```

Remember, point doubling is:

```
x3 = ((3*x1^2 + a)/(2*y1))^2 - 2*x1
y3 = ((3*x1^2 + a)/(2*y1))(x1-x3) - y1
```

Now we move the equations into projective coordinates:

```
x3/(z3^2) = ... = ((3*x1^2 + a * z1^4)^2 - 8 * x1 * y1^2) / (4 * y1^2 * z1^2)
y3/(z3^3) = ... (in the denominators there are y1, z1, z1^2, z1^3)
```

If we choose z3 = 2 * y1 * z1, we obtain formulas for x3 and y3 without any inversions.

Thus, doubling the point (and also addition of two points) in Jacobian coordiantes can be done without modular inversion and is thus much faster as in affine coordinates. Naturally, at the end we have to transform the result back to affine coordinates (which involves modular inversion), but in elliptic curve cryptography the points have to usually be added many times (we have to calculate for example m*Q = Q + ... + Q) and we need to transform to affine coordinates only after last addition.

## Choosing groups and hashing problem

Let E be an elliptic curve over a field K. Every point on E generates (using group operation as above) a cyclic group G. Elliptic curve cryptography is using points in E(K) instead of some F\* for some finite field K, however it needs to be ensured that randomly chosen points and hashed points lie in G, not in all E(K) (see a discussion in [1]).

Cryptographic schemes in F\* for some finite field F operate within a subgroup G of a particular order r, so elements chosen at random and hashed to must have order r or a factor or r. To ensure this, after choosing or hashing to some element x from F\*, to obtain an element from G, x is exponentiated by n/r where n = #F\*. Similarly, in elliptic curve cryptography, P from E(K) is multiplied by n/r where n = #E(K).

If n = #E(K), r is prime that divides n, and r^2 does not divide n, we know there is exactly one subgroup G of E(K) of order r. When a random group element of G is required, we choose a random point of E(K) and multiply it by n/r.When hashing to a point of G, we first hash it to a point in E(K) and then multiply it by n/r.

## Elliptic Curve Discrete Logarithm Problem (ECDLP)

The security of elliptic curve cryptosystems relies on the difficulty of ECDLP:

For two points Q, P from E(F_p) such that Q = k * P, compute k.

The best algorithm for solving ECDLP is Pollard rho method (see ecc_attacks.md).

## Group order

The number of points in E(F_p) is called an order (the order of a finite field is the number of elements it contains, the order of a finite field is always a prime or a power of a prime). The order is denoted as #E(F_p).


Let us choose a point from the curve and check its order:

```
sage: G = F.point((1,2,1))
sage: G.order()
8
```

This means that when multiplying scalars with G there will be 8 different points (so not all points on this curve -> G is not a generator). 

```
sage: [ G * x for x in range(G.order()) ]
[(0 : 1 : 0),
 (1 : 2 : 1),
 (18 : 15 : 1),
 (9 : 12 : 1),
 (7 : 0 : 1),
 (9 : 7 : 1),
 (18 : 4 : 1),
 (1 : 17 : 1)]
```

We can get generators by calling gens():

```
sage: F.gens()
[(3 : 15 : 1)]
sage: G1 = F((3,15,1))
sage: G1.order
24
sage: G1*0
(0 : 1 : 0)
sage: G1*1
(3 : 15 : 1)
sage: G1*2
(5 : 9 : 1)
```

Now, let us compute the inverse of G1:

```
sage: G1 * (-1)
(3 : 4 : 1)
```

It has the same x coordinate. If we sum G1 and its inverse, we shall get infinity (identity):

```
sage: G1 + G1 * (-1)
(0 : 1 : 0)
```

Weierstrass equation has at most two distinct y for each x from F_p, thus #E(F_p) is between 1 and 2*p+1. Hasse's theorem provides some better estimation:

```
(q^(1/2) - 1)^2 <= #E(K) <= (q^(1/2) + 1)^2
```


## Choosing domain parameters

When doing elliptic curve cryptography we have to choose a base point G. This point has to be chosen such that the subgroup generated by it will have a 'sufficiently large' order.

The elliptic curve domain parameters are a sextuple (p, a, b, G, n, h) where: 

 * p specifies the finite field F_p
 * elements a and b specify an elliptic curve E(F_p) defined by the equation E: y^2 = x^3 + a*x + b (mod p)
 * G is a base point on E(F_p)
 * prime n is the order of G (most of the standardized curves have n as a prime number, for example ECDSA cannot be used if n is not prime)
 * h is cofactor, calculated as #E(F_p)/n, where #E(F_p) is number of points on E(F_p)


Domain parameters must be carefully chosen. There exists many elliptic curves or underlying finite fields with special properties which make cryptography over these curves insecure.

See for example Pohlig-Hellman attack in ecc_attacks.md why n (order of G) must be divisible by at least one 'large enough' prime.

### Finding a base point

One way to find a base point is to define some large enough prime n which divides #E(F_p) and search for a point which generates a subgroup of order n. If we choose a random point P on E(F_p) and multiply it with #E(F_p), we get 0:

``` #E(F_p) * P = n * (h*P) = 0 ```

By Lagrange theorem (see [security_definitions.md](https://github.com/miha-stopar/crypto-notes/blob/master/security_definitions.md)) the order of every subgroup H of G divides the order of G. So #E(F_p) is a multiple of order of <P> (subgroup generated by P) and #E(F_p) * P is 0 (because if we multiply P with order of <P> we by definition get 0). 

From the equation above h\*P (or actually the subgroup it generates) is of order n. If we were unlucky and got h*P = 0, we choose another random P and repeat the process.

The value n has to be prime, otherwise we might end up with a point h*P having an order actually lesser than n (one of its divisors).

Elliptic curve cryptograph allows to use key pairs where the lengths of both (public and private) keys are relatively small.


[1] Lynn, Ben. On the implementation of pairing-based cryptosystems. Diss. Stanford University, 2007.

