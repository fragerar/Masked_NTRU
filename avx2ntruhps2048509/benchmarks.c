#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "masking.h"
#include "cpucycles.h"
#include "random.h"
#include "owcpa.h"
#include "params.h"
#include "api.h"

#ifndef ITER
#define ITER 100
#endif

uint64_t start, stop;


void benchmark_masked_poly_map_mod3_to_modq(){
  masked_poly a,b;
  for(int i=0; i < NTRU_N; ++i) for(int j=0; j < N_SHARES; ++j) a.shares[j].coeffs[i] = rand()%3;

  start = cpucycles();
  for(int i=0; i < ITER; i++) masked_poly_map_mod3_to_modq(&b, &a);
  stop = cpucycles();
  printf("Avg speed masked_poly_map_mod3_to_modq: %f cycles.\n", (double)(stop-start)/(ITER));
}

void benchmark_halfmasked_poly_Rq_mul(){
  poly a;
  masked_poly b, c;

  start = cpucycles();
  for(int i=0; i < ITER; i++) halfmasked_poly_Rq_mul(&c, &a, &b); 
  stop = cpucycles();
  printf("Avg speed halfmasked_poly_Rq_mul: %f cycles.\n", (double)(stop-start)/(ITER));
}

void benchmark_masked_poly_mod3_reduce(){
  masked_poly a,b;
  for(int i=0; i < NTRU_N; ++i) for(int j=0; j < N_SHARES; ++j) a.shares[j].coeffs[i] = rand()%NTRU_Q;

  start = cpucycles();
  for(int i=0; i < ITER; i++) masked_poly_mod3_reduce(&b, &a);
  stop = cpucycles();
  printf("Avg speed masked_poly_mod3_reduce: %f cycles.\n", (double)(stop-start)/(ITER));
}

void benchmark_sec_S3_mul(){
  masked_poly a, b, c;
  for(int i=0; i < NTRU_N; ++i) for(int j=0; j < N_SHARES; ++j) a.shares[j].coeffs[i] = rand()%3;
  for(int i=0; i < NTRU_N; ++i) for(int j=0; j < N_SHARES; ++j) b.shares[j].coeffs[i] = rand()%3;

  for(int j=0; j < N_SHARES; ++j) a.shares[j].coeffs[NTRU_N-1] = 0;
  for(int j=0; j < N_SHARES; ++j) b.shares[j].coeffs[NTRU_N-1] = 0;

  start = cpucycles();
  for(int i=0; i < ITER; i++) sec_S3_mul(&c, &a, &b);
  stop = cpucycles();
  printf("Avg speed sec_S3_mul: %f cycles.\n", (double)(stop-start)/(ITER));
}



void benchmark_poly_masked_check_message_space(){
  masked_poly mp;
  poly p;
  Masked b1;
  int cpt=0, index;

  for(int i=0; i < NTRU_N; ++i) {
    p.coeffs[i] = 0;
  }


  while (cpt < NTRU_WEIGHT/2){
    index = rand16()%NTRU_N;
    if (p.coeffs[index] == 0){
      p.coeffs[index] = 1;
      cpt++;
    }
  }
  cpt = 0;
  while (cpt < NTRU_WEIGHT/2){
    index = rand16()%NTRU_N;
    if (p.coeffs[index] == 0){
      p.coeffs[index] = 2;
      cpt++;
    }

  }

  mask_arith_poly(&mp, &p, 3);

  start = cpucycles();
  for(int i=0; i < ITER; i++) poly_masked_check_message_space(&b1, &mp);
  stop = cpucycles();
  printf("Avg speed poly_masked_check_message_space: %f cycles.\n", (double)(stop-start)/(ITER));


}

void benchmark_poly_masked_ternary_check(){
  masked_poly mp;
  poly p;
  Masked t, b;
  for(int i=0; i < NTRU_N; ++i) {
    p.coeffs[i] = rand16()%3;
    if (p.coeffs[i] == 2) p.coeffs[i] = NTRU_Q - 1;
  }
  mask_arith_poly(&mp, &p, NTRU_Q);

  start = cpucycles();
  for(int i=0; i < ITER; i++) poly_masked_ternary_check(&b, &mp);
  stop = cpucycles();
  printf("Avg speed poly_masked_ternary_check: %f cycles.\n", (double)(stop-start)/(ITER));
  
 


}

void benchmark_masked_poly_lift(){
  
  masked_poly p1, p2;
  for(int i=0; i < NTRU_N; ++i) for(int j=0; j < N_SHARES; ++j) p1.shares[j].coeffs[i] = rand()%3;

  start = cpucycles();
  for(int i=0; i < ITER; i++) masked_poly_lift(&p2, &p1);
  stop = cpucycles();
  printf("Avg speed masked_poly_lift: %f cycles.\n", (double)(stop-start)/(ITER));
}

void benchmark_masked_pack_S3(){
  uint8_t packed[NTRU_PACK_TRINARY_BYTES*(N_SHARES)];
  masked_poly p;
  for(int i=0; i < NTRU_N; ++i) for(int j=0; j < N_SHARES; ++j) p.shares[j].coeffs[i] = rand()%3;

  for(int j=0; j < N_SHARES; ++j) p.shares[j].coeffs[NTRU_N-1] = 0;
  
  start = cpucycles();
  for(int i=0; i < ITER; i++) masked_pack_S3(packed, &p);
  stop = cpucycles();
  printf("Avg speed masked_pack_S3: %f cycles.\n", (double)(stop-start)/(ITER));

}

void benchmark_masked_pack_s3_from_Rq(){
  uint8_t packed[NTRU_PACK_TRINARY_BYTES*(N_SHARES)];
  masked_poly p;
  for(int i=0; i < NTRU_N; ++i) for(int j=0; j < N_SHARES; ++j) p.shares[j].coeffs[i] = rand()%NTRU_Q;

  for(int j=0; j < N_SHARES; ++j) p.shares[j].coeffs[NTRU_N-1] = 0;
  
  start = cpucycles();
  for(int i=0; i < ITER; i++) pack_s3_from_Rq(packed, &p);
  stop = cpucycles();
  printf("Avg speed masked_pack_s3_from_Rq: %f cycles.\n", (double)(stop-start)/(ITER));
}


void benchmark_inverse_s3(){

  poly x;
  masked_poly my, mx;

  for(int i=0; i < NTRU_N; ++i) x.coeffs[i] = rand16()%3;
  mask_arith_poly(&mx, &x, 3);


  start = cpucycles();
  for(int i=0; i < ITER; i++) masked_inverse_s3(&my, &mx);
  stop = cpucycles();
  printf("Avg speed masked_inverse_s3: %f cycles.\n", (double)(stop-start)/(ITER));
}

void benchmark_inverse_sq(){

  poly x;
  masked_poly my, mx;

  for(int i=0; i < NTRU_N; ++i) x.coeffs[i] = rand16()%NTRU_Q;
  mask_arith_poly(&mx, &x, NTRU_Q);

  masked_inverse_Sq(&my, &mx);


  start = cpucycles();
  for(int i=0; i < ITER; i++) masked_inverse_Sq(&my, &mx);
  stop = cpucycles();
  printf("Avg speed masked_inverse_Sq: %f cycles.\n", (double)(stop-start)/(ITER));
}





void profiling_decaps(){

  benchmark_sec_S3_mul();
  benchmark_halfmasked_poly_Rq_mul();
  benchmark_masked_poly_mod3_reduce();
  benchmark_poly_masked_ternary_check();
  benchmark_masked_pack_S3();
  benchmark_masked_pack_s3_from_Rq();
  benchmark_poly_masked_check_message_space();
  benchmark_masked_poly_map_mod3_to_modq();
  benchmark_masked_poly_lift();
}

void benchmark_unmasked_owcpa(){
  uint8_t pk[NTRU_PUBLICKEYBYTES]; 
  uint8_t sk[NTRU_SECRETKEYBYTES]; 
  uint8_t c[NTRU_CIPHERTEXTBYTES]; 
  uint8_t k[NTRU_SHAREDKEYBYTES];
  crypto_kem_keypair(pk, sk);
  crypto_kem_enc(c, k, pk);
  crypto_kem_dec(k, c, sk);

}

void benchmark_unmasked(){

  uint8_t pk[NTRU_PUBLICKEYBYTES]; 
  uint8_t sk[NTRU_SECRETKEYBYTES]; 
  uint8_t c[NTRU_CIPHERTEXTBYTES]; 
  uint8_t k[NTRU_SHAREDKEYBYTES];
  crypto_kem_keypair(pk, sk);
  crypto_kem_enc(c, k, pk);

  
  start = cpucycles();
  for(int i=0; i < ITER; i++) crypto_kem_dec(k, c, sk);
  stop = cpucycles();
  printf("Avg speed unmasked decaps: %f cycles.\n", (double)(stop-start)/(ITER));

}

void benchmark_decaps(){

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
  
  masked_crypto_kem_enc(c, k1, pk);
 




  start = cpucycles();
  for(int i=0; i < ITER; i++) masked_crypto_kem_dec(k2, c, &mf, &mf_p, &invh, msk_seed);
  stop = cpucycles();
  printf("Avg speed decaps: %f cycles.\n", (double)(stop-start)/ITER);
}


void benchmark_keygen(){
  masked_poly mf, mf_p;
  poly h_q;
  uint8_t pk[NTRU_OWCPA_PUBLICKEYBYTES];
  uint8_t msk_seed[N_SHARES*32];


  start = cpucycles();
  for(int i=0; i < 10; i++) masked_owcpa_keypair(&mf, &mf_p, &h_q, pk, msk_seed);  
  stop = cpucycles();
  printf("Avg speed keypair: %f cycles.\n", (double)(stop-start)/10);

}


void benchmark_unmasked_keygen(){


  uint8_t pk[NTRU_PUBLICKEYBYTES]; 
  uint8_t sk[NTRU_SECRETKEYBYTES]; 

  start = cpucycles();
  for(int i=0; i < ITER; i++) crypto_kem_keypair(pk, sk);  
  stop = cpucycles();
  printf("Avg speed unmasked keypair: %f cycles.\n", (double)(stop-start)/ITER);

}


int main(){
  printf("Order: %i, rng mode: %i\n", MASKING_ORDER, RNG_MODE);
  //benchmark_inverse_s3();
  //benchmark_inverse_sq();
  benchmark_unmasked();
  benchmark_unmasked_keygen();
  benchmark_decaps();
  benchmark_keygen();
  profiling_decaps();
  printf("\n");
  return 0;
}
