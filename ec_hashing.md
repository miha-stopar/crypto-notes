# EC hashing

## Tonelli-Shanks algorithm

Tonelli-Shanks algorithm is used to solve r^2 = n (mod p), where p is prime.

Note that if n has a square root modulo p, it holds (Euler criterion): 

```
n^((p-1)/2) = 1 (mod p). 
```

In this case we can write (we omit modulo):

```
n = n^((p-1)/2) * n = n^((p+1)/2)
```

Thus x = n^((p+1)/4)) is a square root of n. However, this is true only if p+1 is divisible by 4.

Note that we have equation of the form:

```
x^2 = t * n.
```

When p+1 is not divisible by 4, we use Tonelli-Shanks to modify this equation to:

```
x1^2 = t1 * n 
```

where the order of t1 is smaller than the order of t (and thus gradually converge to t1 with order 0, which is 1).

First, we modify the initial equation to (e is the biggest such that p-1 = Q * 2^e):

```
x^2 = (n^(((p-1)/2^e + 1)/2)^2 = n^((p-1)/2^e) * n
```

Note that `(p-1)/2^e + 1` is even, so we can divide it by 2.
If n^((p-1)/2^e) = 1, we are done. Otherwise we find the lowest i such that:

```
(n^((p-1)/2^e))^2^i = 1
```

We find z which does not have a square root modulo p (quadratic nonresidue). Element z^((p-1)/2^i) has order 2^i (it can't be less because z^((p-1)/2 = -1 by Euler criterion).

We show that z^((p-1)/2^i) is quadratic residue:

```
z^((p-1)/2^i) = (z^((p-1)/2^e)^2^(e-i-1))^2
```

We multiply the equation by z (let's say z1^2 = z).

```
x^2 = t * n
x^2 * z1^2 = t * z1^2 * n
```

We have t1 = t * z1^2 with order smaller than order of t.

## Generalized Tonelli-Shanks algorithm

Let's compute a^(1/l) mod p.

```
p-1 = l^r * u where l is prime and it does not divide u
```

Let's compute x, y such that: `y*l - u*x = 1`. So:

```
a^(y*l - u*x) = a
(a^y)^l = a^(u*x) * a
```

If a^(u*x) = 1, we have a solution: a^y. If not, we multiply the equation by some l-power element b (b = b1^l) which makes the order of a^(u*x) * b smaller than the order of a^(u*x).

Let's say the order of a^(u*x) is l^w (we know it is a power of l). Now we choose a generator g of a l-Sylow group. It holds: g^(l^r) = 1. We find v such that: g^(l^v) has the same order as a^(u*x).

If we apply homomorphism f: x -> x^(l^(w-1)) to a^(u*x) and g^(l^v), we get two l-roots. That means both are from the same cyclic group and we can find z such that:

```
f(a^(u*x)) = f(g^(l^v))^z
f(a^(u*x)) = f(g^(z * l^v))
```

Now the element a^(u*x) * g^(-z * l^v) is of order l^(w-1) because f(a^(u*x) * g^(-z * l^v)) = 1. Thus we reduced the order and we repeat the process until the order is 0. Note that g^(z * l^v) is a power of l, so the left side of the equation is still a power of l.

In a special case when r = 1 (p-1 = l^1 * u), we can write for the l-th power a = x^l:

```
a^((p-1)/l) = 1
a^((p-1+l)/l) = a
a^((l*u+l)/l) = a
a^(u+1) = a
```

If u+1 is divisible with l, we can easily obtain l-th root of a: a^((u+1)/l).


