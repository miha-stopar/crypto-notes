#

Legend:

```
isomorphic: ~= 
```

Projective space over a field K is a set of one-dimensional subspaces of the vector space K^(n+1) and is denoted as P^n. Sometimes, the projective space (of one-dimensional subspaces) of a vector space V is denoted by P(V).

Example: Linear spaces

An inclusion of vector spaces W ~= K^(k+1) -> V ~= K^(n+1) induces a map P(W) -> P(V). The image lambda of such a map is called a linear subspace of dimension k in P(V). In case k = n-1, we call it a hyperplane.

If we consider a span of two linear subspaces lambda = P(W), lambda1 = P(W1): this is the subspace associated to the sum W + W1 (the smallest linear subspace containing both lambda and lambda1), denoted by lambda,lambda1.

The dimension of the subspace lambda,lambda1 is at most the sum of the dimensions plus one, with equality holding iff lambda and lambda1 are disjoint. In general:

```
dim(lambda,lambda1) = dim(lambda) + dim(lambda1) - dim(lambda ∩ lambda1)
```

where we take the dimension of the empty set as a linear subspace to be -1. One of the basic properties of projective space P^n: whenever k+l >= n, any two linear subspaces lambda, lambda1 of dimension k and l in P^n will intersect in a linear subspace of dimension at least k+l-n.

Every hyperplane H in P^n can be expressed as the kernel of a non-trivial linear form, a K-linear map:

```
phi: K^(n+1) -> K
x = (x_0,...,x_n) -> u_0 * x_0 + ... + u_n * x_n
```

The set of all K-linear forms on K^(n+1) yields the dual space (K^(n+1))\*.

