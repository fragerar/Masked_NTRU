#include "masking.h"
#include "poly.h"
#include "owcpa.h"
#include "sample.h"
#include "api.h"
#include "cmov.h"
#include "crypto_hash_sha3256.h"
#include "kem.h"
#include "rng.h"
#include "params.h"

int masked_crypto_kem_enc(unsigned char *c, unsigned char *k, const unsigned char *pk)
{
  masked_poly mr, mm;
  unsigned char mrm[NTRU_OWCPA_MSGBYTES*(N_SHARES)];
  uint8_t packed_r[NTRU_PACK_TRINARY_BYTES*(N_SHARES)], packed_m[NTRU_PACK_TRINARY_BYTES*(N_SHARES)];
  


  masked_ternary(&mr);
  #ifdef NTRU_HPS
  masked_fixed_type(&mm, 3);
  #endif

  #ifdef NTRU_HRSS
  masked_ternary(&mm);
  #endif

  masked_pack_S3(packed_r, &mr);
  masked_pack_S3(packed_m, &mm);


  for(int i=0; i < NTRU_PACK_TRINARY_BYTES; ++i){
    for(int j=0; j < N_SHARES; ++j){
      mrm[i                          + j*2*NTRU_PACK_TRINARY_BYTES] = packed_r[i + j*NTRU_PACK_TRINARY_BYTES];
      mrm[i+ NTRU_PACK_TRINARY_BYTES + j*2*NTRU_PACK_TRINARY_BYTES] = packed_m[i + j*NTRU_PACK_TRINARY_BYTES];
    }
  }

  sha3_256_masked(k, mrm, NTRU_OWCPA_MSGBYTES);

  masked_poly_map_mod3_to_modq(&mr, &mr);
  masked_owcpa_enc(c, &mr, &mm, pk);

  return 0;
}

int masked_crypto_kem_dec(unsigned char *k, const unsigned char *c, const masked_poly* mf, const masked_poly* mf_p, const poly* invh, const uint8_t* msk_seed)
{
  int i, fail;
  unsigned char mrm[NTRU_OWCPA_MSGBYTES*(N_SHARES)];
  unsigned char mbuf[(NTRU_PRFKEYBYTES+NTRU_CIPHERTEXTBYTES)*N_SHARES];

  fail = masked_owcpa_dec(mrm, c, mf, mf_p, invh);
  /* If fail = 0 then c = Enc(h, rm). There is no need to re-encapsulate. */
  /* See comment in owcpa_dec for details.                                */

  sha3_256_masked(k, mrm, NTRU_OWCPA_MSGBYTES);

  /* shake(secret PRF key || input ciphertext) */  
  for(int j=0; j < N_SHARES; ++j){
    for(i=0;i<NTRU_PRFKEYBYTES;i++)
      mbuf[i + j*(NTRU_PRFKEYBYTES+NTRU_CIPHERTEXTBYTES)] = msk_seed[i + j*NTRU_PRFKEYBYTES];
    for(i=0;i<NTRU_CIPHERTEXTBYTES;i++)
      mbuf[NTRU_PRFKEYBYTES + i + j*(NTRU_PRFKEYBYTES+NTRU_CIPHERTEXTBYTES)] = 0;
  }
  for(i=0;i<NTRU_CIPHERTEXTBYTES;i++) mbuf[NTRU_PRFKEYBYTES + i] = c[i];

  sha3_256_masked(mrm, mbuf, NTRU_PRFKEYBYTES+NTRU_CIPHERTEXTBYTES);
  
  cmov(k, mrm, NTRU_SHAREDKEYBYTES*(N_SHARES), (unsigned char) fail);

  return 0;
}