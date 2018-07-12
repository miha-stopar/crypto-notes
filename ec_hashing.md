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





