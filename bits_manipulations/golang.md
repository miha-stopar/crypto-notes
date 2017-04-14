Convert string to ASCII codes:
```
h := []byte("hello")
fmt.Println(h)
[104 101 108 108 111]
```

From bytes to string:

```
fmt.Println(string(h)
hello
```

Convert to an integer (import "math/big"):
```
fmt.Println(new(big.Int).SetBytes(h))
448378203247 //  104*256**4 + 101*256**3 + 108*256**2 + 108*256 + 111
```

Create a large integer:

```
N, _ := new(big.Int).SetString("26959946667150639794667015087019630673557916260026308143510066298881", 10)
```

Back to string:

```
N.String()
```

Or from hex:

```
B, _ := new(big.Int).SetString("b4050a850c04b3abf54132565044b0b7d7bfd8ba270b39432355ffb4", 16)
```

Transform into bytes:

```
fmt.Println(B.Bytes())
[180 5 10 133 12 4 179 171 245 65 50 86 80 68 176 183 215 191 216 186 39 11 57 67 35 85 255 180]

```











