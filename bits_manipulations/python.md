Convert character to ASCII:
```
In [1]: ord('a')
Out[1]: 97
```

Convert ASCII value to character:
```
In [2]: chr(97)
Out[2]: 'a'
```

Convert string to hex:
```
In [3]: "hello".encode("hex")
Out[3]: '68656c6c6f'
```

And back:
```
In [4]: "68656c6c6f".decode("hex")
Out[4]: 'hello'
```

Convert hex to number:
```
In [5]: int('68656c6c6f', 16)
Out[5]: 448378203247
In [6]: num = 448378203247
```

To convert an integer to a binary string of bits:
```
In [7]: bin_num = bin(num)
In [8]: bin_num
Out[8]: '0b110100001100101011011000110110001101111'
```

To see the number of bits (you have to remove 2 0 for initial '0b'):
```
In [10]: len(bin_num)
Out[10]: 41

In [11]: num.bit_length()
Out[11]: 39
```

Convert back to integer:
```
In [12]: int(bin_num, 2)
```

Convert integer to hex:
```
In [13]: hex(num)
Out[13]: '0x68656c6c6f'
```

There is a module named struct which can be of great help when converting between strings and bytes. For example if you have a message m = "12345678" * 8 = "1234567812345678123456781234567812345678123456781234567812345678" and you want to convert it into 16 numbers where each number will represent 4 characters (16 * 4 = 64):

```
In [26]: struct.unpack("<16I", m)
Out[26]: 
(875770417,
 943142453,
 875770417,
 943142453,
 875770417,
 943142453,
 875770417,
 943142453,
 875770417,
 943142453,
 875770417,
 943142453,
 875770417,
 943142453,
 875770417,
 943142453)
```

Note that 875770417 represents "1234" and 943142453 represents "5678". The first argument "<16I" specifies that we want to convert into 16 ("16") unsigned integers ("I") using little-endian ("<") byte order (the most significant bit is the last one). That means we take for example "1234", check the ASCII value for each of the ["1", "2", "3", "4"] which are [49, 50, 51, 52] and calculate  52*256**3 + 51*256**2 + 50*256 + 49 = 875770417.




