#ifndef API_H
#define API_H


#ifdef PARAMS_1
  #define CRYPTO_SECRETKEYBYTES 935
  #define CRYPTO_PUBLICKEYBYTES 699
  #define CRYPTO_CIPHERTEXTBYTES 699
  #define CRYPTO_BYTES 32

  #define CRYPTO_ALGNAME "ntruhps2048509"
#endif

#ifdef PARAMS_2
  #define CRYPTO_SECRETKEYBYTES 1234
  #define CRYPTO_PUBLICKEYBYTES 930
  #define CRYPTO_CIPHERTEXTBYTES 930
  #define CRYPTO_BYTES 32

  #define CRYPTO_ALGNAME "ntruhps2048677"
#endif

#ifdef PARAMS_3
  #define CRYPTO_SECRETKEYBYTES 1450
  #define CRYPTO_PUBLICKEYBYTES 1138
  #define CRYPTO_CIPHERTEXTBYTES 1138
  #define CRYPTO_BYTES 32

  #define CRYPTO_ALGNAME "ntruhrss701"
#endif

#ifdef PARAMS_4
  #define CRYPTO_SECRETKEYBYTES 1590
  #define CRYPTO_PUBLICKEYBYTES 1230
  #define CRYPTO_CIPHERTEXTBYTES 1230
  #define CRYPTO_BYTES 32

  #define CRYPTO_ALGNAME "ntruhps4096821"
#endif


#define crypto_kem_keypair CRYPTO_NAMESPACE(keypair)
int crypto_kem_keypair(unsigned char *pk, unsigned char *sk);

#define crypto_kem_enc CRYPTO_NAMESPACE(enc)
int crypto_kem_enc(unsigned char *c, unsigned char *k, const unsigned char *pk);

#define crypto_kem_dec CRYPTO_NAMESPACE(dec)
int crypto_kem_dec(unsigned char *k, const unsigned char *c, const unsigned char *sk);

#endif
