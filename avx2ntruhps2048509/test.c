#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include "poly.h"
#include "owcpa.h"
#include "masking.h"
#include "random.h"
#include "api.h"
#include "sample.h"
#include "rng.h"


void random_poly(poly *p, int q){
  for(int i=0; i < NTRU_N; ++i) p->coeffs[i] = rand()%q;
}


void test_masked_kem(){
  uint8_t pk[NTRU_PUBLICKEYBYTES]; 
  uint8_t sk[NTRU_SECRETKEYBYTES]; 
  uint8_t c[NTRU_CIPHERTEXTBYTES]; 
  uint8_t k1[NTRU_SHAREDKEYBYTES*(N_SHARES)];
  uint8_t k2[NTRU_SHAREDKEYBYTES*(N_SHARES)];
  uint8_t sk_seed[NTRU_PRFKEYBYTES];
  uint8_t msk_seed[NTRU_PRFKEYBYTES*(N_SHARES)];
  

  poly f, finv3, invh;
  masked_poly mf, mf_p;

  crypto_kem_keypair(pk, sk);

  poly_S3_frombytes(&f, sk);
  poly_S3_frombytes(&finv3, sk+NTRU_PACK_TRINARY_BYTES);
  poly_Sq_frombytes(&invh, sk+2*NTRU_PACK_TRINARY_BYTES);
  mask_arith_poly(&mf, &f, 3);
  mask_arith_poly(&mf_p, &finv3, 3);
  for(int i=0;i<NTRU_PRFKEYBYTES;i++) sk_seed[i] = sk[i+NTRU_OWCPA_SECRETKEYBYTES];
  mask_bitstring(msk_seed, sk_seed, NTRU_PRFKEYBYTES);
  
/*   masked_crypto_kem_enc(c, k1, pk);
  printf("k1 : "); print_masked_bs(k1, NTRU_SHAREDKEYBYTES); */

  uint8_t k[NTRU_SHAREDKEYBYTES];
  crypto_kem_enc(c, k, pk);
  printf("k1 : "); for(int i=0; i < NTRU_SHAREDKEYBYTES; ++i) printf("%X ", k[i]); printf("\n");


  masked_crypto_kem_dec(k2, c, &mf, &mf_p, &invh, msk_seed);
  printf("k2 : "); print_masked_bs(k2, NTRU_SHAREDKEYBYTES);


}


void test_masked_keypair(){

  masked_poly mf, mf_p;
  poly h_q;
  uint8_t pk[NTRU_OWCPA_PUBLICKEYBYTES];
  uint8_t msk_seed[N_SHARES*32];
  uint8_t c[NTRU_CIPHERTEXTBYTES]; 
  uint8_t k2[NTRU_SHAREDKEYBYTES*(N_SHARES)];

  masked_owcpa_keypair(&mf, &mf_p, &h_q, pk, msk_seed);


  uint8_t k[NTRU_SHAREDKEYBYTES];
  crypto_kem_enc(c, k, pk);
  printf("k1 : "); for(int i=0; i < NTRU_SHAREDKEYBYTES; ++i) printf("%X ", k[i]); printf("\n");
  

  masked_crypto_kem_dec(k2, c, &mf, &mf_p, &h_q, msk_seed);

  printf("k2 : "); print_masked_bs(k2, NTRU_SHAREDKEYBYTES);

}

  
int main(){
  //srand(1);
  srand(time(0));

  unsigned char entropy[48], pstring[48];
  for(int i=0; i < 48; ++i) entropy[i] = rand();
  randombytes_init(entropy, pstring, 128);

  printf("Hello world\n");

  test_masked_keypair(); printf("\n\n");
  test_masked_kem(); printf("\n\n");

  

  
  printf("=====\n");



  return 0;
}