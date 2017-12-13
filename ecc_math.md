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
4 * y^2 + 4 * a_1 * x * y + 4 * a_3 * y = 4 * x^3 + 4 * a_2 * x^2 + a_4 * x + a_6
(2 * y + a_1 * x + a_3)^2 - a_1^2 * x^2 - 2 * a_1 * a_3 * x - a_3^2 = 4 * x^3 + 4 * a_2 * x^2 + 4 * a_4 * x + 4 * a_6
(2 * y + a_1 * x + a_3)^2 = 4 * x^3 + (4 * a_2 + a_1^2) * x^2 + (4 * a_4 + 2 * a_1 * a_2) * x + (4 * a_6 + a_3^2)
```

We get:

```
y^2 = 4 * x^3 + b_2*x^2 + 2 * b_4*x + b_6
```

where:

```
b_2 = a_1^2 + 4 * a_2
b_4 = a_1*a_3 + 2 * a_4
b_6 = a_3^2 + 4 * a_6
```

We define the following quantities:

```
b_8 = a_1^2 * a_6 + 4 * a_2 * a_6 - a_1 * a_3 * a_4 + a_2 * a_3^2 - a_4^2
c_4 = b_2^2 - 24 * b_4
c_6 = -b_2^3 + 36 * b_2 * b_4 - 216 * b_6
delta = -b_2^2 * b_8 - 8 * b_4^3 - 27 * b_6^2 + 9 * b_2 * b_4 * b_6
j = c_4^3/delta
omega = dx / (2 * y + a_1 * x + a_3) = dy / (3 * x^2 + 2 * a_2 * x + a_4 - a_1 * y)
```

By definition delta is the discriminant, j is the j-invariant, w is the invariant differential. It turns out if delta = 0, the curve is singular.

It holds:

```
4 * b_8 = b_2 * b_6 - b_4^2
1728 * delta = c_4^3 - c_6^2
```

If char(K~) != 2, 3, then with additional substition:

```
(x,y) -> ((x - 3 * b_2)/36, y/108)
```

we eliminate x^2 term and get:

```
E: y^2 = x^3 - 27 * c_4 * x - 54 * c_6
```

Let's have a look at the singular points of E. We have f(x, y):

```
f(x, y) = y^2 + a_1*x*y + a_3*y - x^3 - a_2*x^2 - a_4*x - a_6
```

At singular point P = (x0, y0):

```
(df/dx)(P) = (df/dy)(P) = 0
```

Taylor expansion at P:

```
f(x, y) = f(x0, y0) + (df/dxdx)(P) * ((x-x0)^2)/2 + (df/dxdy)(P) * (x-x0)*(y-y0)/2 + (df/dydy)(P) * ((y-y0)^2)/2 + (df/dxdxdx)(P) * ((x-x0)^3)/6
```

All other derivatives are 0 ((df/dxdxdy), (df/dxdydy), (df/dydydy)), also f(x0, y0) = 0 and (df/dxdxdx)(P) * ((x-x0)^3)/6 = (x-x0)^3:

```
f(x, y) = (df/dxdx)(P) * ((x-x0)^2)/2 + (df/dxdy)(P) * (x-x0)*(y-y0)/2 + (df/dydy)(P) * ((y-y0)^2)/2 + (x-x0)^3
```

The tangent is:

```
(df/dxdx)(P) * (x-x0)^2 + (df/dxdy)(P) * (x-x0)*(y-y0) + (df/dydy)(P) * (y-y0)^2 = 0
```

We can get alpha, beta such that the equation above becomes:

```
((y-y0) - alpha * (x-x0)) * ((y-y0) - beta * (x-x)) = 0
```

If alpha = beta, we have one tangent, otherwise two. When alpha = beta, we say P is a cusp, otherwise a node.


When char(K) is not 2 or 3, the elliptic curve has a Weierstrass equations of the form:

```
y^2 = x^3 + a * x + b
```

and the quantities:

```
delta = -16 * (4 * a^3 + 27 * b^2)
j = -1728 * (4 * a)^3 / delta
```

For the curve given by a Weierstrass equation it holds (Silverman [1], Proposition 1.4):

 * It is nonsingular iff delta != 0 
 * It has a node iff delta = 0 and c_4 != 0
 * It has a node iff delta = c_4 = 0

Note that (in second and third point) there is only one singular point.

For a proof we need to realize that delta and c_4 stay the same when the transformation x1 = x - x0, y1 = y - y0 is applied; that means we can assume that the singular point is at (0, 0) ...








[1] Silverman, Joseph H. The arithmetic of elliptic curves. Vol. 106. Springer Science & Business Media, 2009.

