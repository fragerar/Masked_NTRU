#ifndef MASKED_KEM
#define MASKED_KEM

#include "params.h"
#include "poly.h"
#include "masking.h"

int masked_crypto_kem_enc(unsigned char *c, unsigned char *k, const unsigned char *pk);
int masked_crypto_kem_dec(unsigned char *k, const unsigned char *c, const masked_poly* mf, const masked_poly* mf_p, const poly* invh, const uint8_t* msk_seed);

#endif