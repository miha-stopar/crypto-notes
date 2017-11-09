
# Example of Weil pairing

Let's check elliptic curve E: y^2 = x^3 + x over F_59. Weil pairing on this curve is discussed in [1].

```
sage: E = EllipticCurve(GF(59), [1,0])
sage: E
Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 59
sage: E.cardinality()
60
```

We can see that it has 60 points.



[1] Lynn, Ben. On the implementation of pairing-based cryptosystems. Diss. Stanford University, 2007.


