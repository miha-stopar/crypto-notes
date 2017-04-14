# OpenSSL

## Basic commands

Encrypt "hello world" using AES-256-CTR and some dummy key and iv:

```
echo -n "hello world" | openssl enc -aes-256-ctr -K 0000000000000000000000000000000000000000000000000000000000000000 -iv 00000000000000000000000000000000
```

If you want to also transform into base64:

```
echo -n "hello world" | openssl enc -aes-256-ctr -K 0000000000000000000000000000000000000000000000000000000000000000 -iv 00000000000000000000000000000000 | base 64
```


## How to run HTTPS server using s_server

First, create a certificate and key:
```
openssl req -x509 -newkey rsa:2048 -keyout key.pem -out cert.pem -days 365 -nodes
```

Start server (apps/s_server.c):
```
openssl s_server -key key.pem -cert cert.pem -accept 44330 -www
```

If you locally made modifications of openssl, make sure you are running the right openssl, for example (after execution of './config; make; make test; sudo make install' in the openssl folder where you are making changes): 
```
/usr/local/ssl/bin/openssl s_server -key key.pem -cert cert.pem -accept 44330 -www
```

Now we can access it from the browser on: https://localhost:44330/

Or using s_client (apps/s_client.c):
```
openssl s_client -connect localhost:44330
```

todo - DROWN:
```
openssl s_client -connect localhost:443 -ssl2
```

For more log messages:
```
openssl s_client -connect localhost:44330 -state -debug
```

There is msg_cb callback (apps/s_cb.c) which handles different messages (Handshake, ChangeCipherSpec...).

In s_client.c there is:

ctx = SSL_CTX_new(meth);

con = SSL_new(ctx);

con is then pushed to print_stuff which output info about chosen cipher ...

Inside ssl/record/rec_layer_s3.c records are created (also countermeasure against known-IV weakness in CBC)

ssl/statem/statem.c implements SSL/TLS/DTLS state machines. There are two primary state machines:
 * message flow state machine
 * handshake state machine

Inside ssl/statem/statem_clnt.c there is "Construct a message to be sent from the client to the server"

messages from server to client are constructed in ssl/statem/statem_srvr.c (Klima-Pokorny-Rosa mentioned)

You can also see in statem_srvr.c how countermeasures for Bleichenbacher attack (and Klima-Pokorny-Rosa extension) are implemented. The crucial thing is obviously not to expose what kind of error happened (if any), but also to ensure that the execution time does not differ for different errors that might happen (to avoid side-channel attacks).

For example random premaster secret (rand_premaster_secret) is generated before premaster_secret is decryted (where the error about message not being PKCS#1 compliant might appear and exposed as a decryption oracle) - ensuring everything is the same no matter whether the error occurs or not. Later on the following is executed:

``` 
for (j = 0; j < sizeof(rand_premaster_secret); j++) {
    rsa_decrypt[j] =
        constant_time_select_8(decrypt_good, rsa_decrypt[j],
                               rand_premaster_secret[j]);
}
``` 

This means in the case of non-compliant message, rand_master_secret is copied over rsa_decrypt, but even if the message was compliant, the copying will happen anyway (rsa_decrypt over rsa_decrypt) and will take the same amount of time due to constant_time_select_8 function.



## How to run tests

In root folder for example to run Diffie-Hellman test: 
```
make TESTS="test_dh" test
```

This will run 15-test_dh.t from openssl/test/recipes. Inside it there is a call:
```
simple_test("test_dh", "dhtest", "dh");
```



# scrypt

scrypt is a password-based key derivation function. It was designed not to be only deliberately slow (see [crypto_constructs.md](https://github.com/miha-stopar/crypto-notes/blob/master/crypto_constructs.md) about password-based KDFs), but also to require a large amount of memory which makes parallelization harder. It is believed to be more secure than alternative solutions such as PBKDF2 or bcrypt.

# NaCl

NaCl is a crypto library which aims to be so easy that it prevents crypto disasters due to improper use.

Let us say that Alice wants to send a message m to Bob. 
Alice has public key pk1 and secret key sk1, while Bob has public key pk2 and secret key sk2.
Alice computes the ciphertext as:

c = crypto_box(m, n, pk2, sk1)

Parameter n is 24-byte nonce, c is authenticated ciphertext which is 16 bytes longes than m.

Bob verifies the authentication and recovers the message using function crypto_box_open.

m = crypto_box_open(c, n, pk1, sk2)

A pair of a public and a private key is generated using function crypto_box_keypair.

Function crypto_sign is used to sign a message, crypto_sign_open to verify the signature and recover the message. Function crypto_sign_keypair is used to generate a keypair for signature scheme.

NaCl uses uses /dev/urandom as random number generator.

# SJCL

When calculating hash, SJCL can take bitArray or String.

> b=sjcl.codec.utf8String.toBits('hi mom')
[ 1751720045, 17594055458816 ]

The following two calls give the same result:

> sjcl.hash.sha256.hash(b)
[ -185578762,
  -1473069864,
  -1501209099,
  -802935654,
  1911151127,
  -85250765,
  -1291024125,
  -239752892 ]

> sjcl.hash.sha256.hash('hi mom')
[ -185578762,
  -1473069864,
  -1501209099,
  -802935654,
  1911151127,
  -85250765,
  -1291024125,
  -239752892 ]

sjcl mostly operates on arrays of 4 byte words, as you can see above. The words 'fooo' and 'fooofooo' are for example presented as:
> sjcl.codec.utf8String.toBits('fooo')
[ 1718579055 ]
> sjcl.codec.utf8String.toBits('fooofooo')
[ 1718579055, 1718579055 ]

Obviosly, 1718579055 = 111 + 111*256 + 111*256*256 + 102*256*256*256 where:
> 'f'.charCodeAt(0)
102
> 'o'.charCodeAt(0)
111

The problem appears when you have a string of less then four characters. You might assume that for example 'foo' is represented as 111 + 111*256 + 102*256*256 or 111*256 + 111*256*256 + 102*256*256*256. But if for example the first case would be used, the string 'foo' would have the same number as the string '\u0000foo', where:
> a=String.fromCharCode(0)+'foo'
'\u0000foo'
> sjcl.codec.utf8String.toBits(a)
[ 6713199 ]

Thus, sjcl does the trick and shifts the number representation of the word (of less then four characters) to the left (by 8, 16 or 24 bits if the string is of length 24, 16 or 8 bits respectively), and then it adds len * 0x10000000000 (where len is the bit length of the string) to the number, for example:
> tmp = 6713199
6713199
> (tmp << 8) + 24 * 0x10000000000
26389997645568

And this is the same as:
> sjcl.codec.utf8String.toBits('foo')
[ 26389997645568 ]

The value 0x10000000000 is the same as 2**40 (or 16**10).

The point is that the value 26389997645568 still needs to return 6713199 when '>> 8' is applied. So the same as would return 6713199 << 8 (= 1718578944). And it does so, because bitwise operation converts the number to a 32-bit signed integer and discards all higher-placed bits. So for example:

> (Math.pow(2,40) + 10) >> 32
10

If you from some reason have to hash a concatenation of two strings, you can obviously do:

> sjcl.hash.sha256.hash('foo' + 'bar')

However, that is not the same as: 

> f = sjcl.codec.utf8String.toBits('foo')
> b = sjcl.codec.utf8String.toBits('bar')
> sjcl.hash.sha256.hash(f.concat(b))

Using the method toBits you can properly (in sjcl context) concatenate two arrays of words:

c = sjcl.bitArray.concat(f, b)
sjcl.hash.sha256.hash(c) is the same as sjcl.hash.sha256.hash('foobar')


http://www.2ality.com/2012/07/large-integers.html

http://kipirvine.com/asm/workbook/floating_tut.htm
http://floating-point-gui.de/languages/javascript/

var x = 999999999999999;   // x will be 999999999999999
var y = 9999999999999999;  // y will be 10000000000000000



