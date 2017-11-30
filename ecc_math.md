# Arithmetic of elliptic curves

Notes mostly taken from Silverman [1].


The following notation will be used:

```
k~: an algebraic closure of field k
```

## Weierstrass equations

An elliptic curve has an equation of the form:

```
Y^2 * Z + a_1 * X * Y * Z + a_3 * Y * Z^2 = X^3 + a_2 * X^2 * Z + a_4 * X * Z^2 + a_6 * Z^3
```

Here O = [0, 1, 0] is the base point and a_1, a_2, a_3, a_4, a_6 are from K~. Using non-homogeneous coordinates (x = X/Z, y = Y/Z), we get:

```
y^2 + a_1*x*y + a_3*y = x^3 + a_2*x^2 + a_4*x + a_6

```

where we need to remember that there is an extra point O = [0, 1, 0] out at infinity. If a_1, a_2, a_3, a_4, a_6 are from K, we say that E is defined over K.

If the characteristic of K is not 2, we can complete the square:

```
(y + 1/2 * (a_1 * x + a_3))^2 = ...
```

We get:

```
y^2 = x^3 + b_2*x^2 + b_4*x + b_6
```

where:

```
b_2 = 1/4 * a_1^2 + a_2
b_4 = 1/2 * a_1*a_3 + a_4
b_6 = 1/4 * a_3^2 + a_6
```






[1] Silverman, Joseph H. The arithmetic of elliptic curves. Vol. 106. Springer Science & Business Media, 2009.

