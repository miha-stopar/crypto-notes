Convert character to ASCII:
"h".ord
=> 104

Convert ASCII value to character:
104.chr
=> "h"

Convert string to hex:
"hello".unpack("H*")
=> ["68656c6c6f"]

Convert hex to string:
["68656c6c6f"].pack("H*")
=> "hello"

Convert hex to number:
"68656c6c6f".hex
=> 448378203247

So, to convert string (not in hex) to a number presentation:
"hello".unpack("H*")[0].hex
=> 448378203247

The number value can be obviously calculated also:
"o".ord + "l".ord * 256 + "l".ord * 256**2 + "e".ord * 256**3 + "h".ord * 256**4
=> 448378203247

Convert hex string to a number:
irb(main):001:0> a = 18.chr * 3 + 19.chr * 2
=> "\x12\x12\x12\x13\x13"
irb(main):002:0> a.unpack("H*")
=> ["1212121313"]
irb(main):003:0> a.unpack("H*")[0].to_i(16)
=> 77612585747



The following huge number can be for example used as p in Diffie-Hellman.

p = 'ffffffffffffffffc90fdaa22168c234c4c6628b80dc1cd129024e088a67cc74020bbea63b139b22514a08798e3404ddef9519b3cd3a431b302b0a6df25f14374fe1356d6d51c245e485b576625e7ec6f44c42e9a637ed6b0bff5cb6f406b7edee386bfb5a899fa5ae9f24117c4b1fe649286651ece45b3dc2007cb8a163bf0598da48361c55d39a69163fa8fd24cf5f83655d23dca3ad961c62f356208552bb9ed529077096966d670c354e4abc9804f1746c08ca237327ffffffffffffffff'

To convert it into number:
b = p.to_i(16)
puts b

This will print a 384 digits number:
2410312426921032588552076022197566074856950548502459942654116941958108831682612228890093858261341614673227141477904012196503648957050582631942730706805009223062734745341073406696246014589361659774041027169249453200378729434170325843778659198143763193776859869524088940195577346119843545301547043747207749969763750084308926339295559968882457872412993810129130294592999947926365264059284647209730384947211681434464714438488520940127459844288859336526896320919633919

To convert it back into hex:
b.to_s(16)

Number b is too big to calculate for example q**b % p, thus for example openssl could be included, then to_bn method could be used (to convert to OpenSSL::BN) and finally mod_exp:
require 'openssl'
q.to_bn.mod_exp(b, p)





s = [p].pack("H*") # string (obviously unreadable)


