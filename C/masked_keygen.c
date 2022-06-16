#include "poly.h"
#include "masking.h"
#include "random.h"


#if NTRU_N == 509
#include "masked_inverse_509.c"

#elif NTRU_N == 677
#include "masked_inverse_677.c"

#elif NTRU_N == 701
#include "masked_inverse_701.c"

#elif NTRU_N == 821
#include "masked_inverse_821.c"

#endif



void masked_inverse_Sq(masked_poly* y, masked_poly* a){
  masked_poly p1, p2, p3;
  int t=1;
  for(int j=0; j < N_SHARES; ++j) for(int i=0; i < NTRU_N; ++i) (p1.shares[j]).coeffs[i] = ((a->shares[j]).coeffs[i])&1;
  masked_inverse_s2(&p2, &p1);
  masked_poly_map_mod2_to_modq(&p1, &p2);


  while (t <= NTRU_LOGQ){
      poly_refresh(&p1, NTRU_Q);
      sec_Sq_mul(&p3, &p1, a);
      for(int j=0; j < N_SHARES; ++j) for(int i=0; i < NTRU_N; ++i) (p3.shares[j]).coeffs[i] = (NTRU_Q - (p3.shares[j]).coeffs[i])%NTRU_Q;
      p3.shares[0].coeffs[0] = (p3.shares[0].coeffs[0] + 2)%NTRU_Q; 
      sec_Sq_mul(&p2, &p3, &p1);
      for(int j=0; j < N_SHARES; ++j) for(int i=0; i < NTRU_N; ++i) (p1.shares[j]).coeffs[i] = ((p2.shares[j]).coeffs[i])%NTRU_Q;
      t = 2*t;
  }
  for(int j=0; j < N_SHARES; ++j) for(int i=0; i < NTRU_N; ++i) (y->shares[j]).coeffs[i] = (p1.shares[j]).coeffs[i];
}


void masked_owcpa_keypair(masked_poly* mf, masked_poly* mf_p, poly* h_q, uint8_t* pk, uint8_t* msk_seed){
  poly h;
  masked_poly mg, mf_q, x, y;

#ifdef NTRU_HPS
  masked_ternary(mf);
  masked_fixed_type(&y, 3);

  masked_poly_map_mod3_to_modq(&mg, &y);
#endif

#ifdef NTRU_HRSS
  masked_ternary_plus(mf);
  masked_ternary_plus(&mg);


  masked_poly_map_mod3_to_modq(&y, &mg);


  for(int j=0; j < N_SHARES; ++j){
    for(int i=NTRU_N-1; i>0; i--)
      mg.shares[j].coeffs[i] = (y.shares[j].coeffs[i-1] + NTRU_Q - y.shares[j].coeffs[i])%NTRU_Q;
    mg.shares[j].coeffs[0] = NTRU_Q-(y.shares[j].coeffs[0]);
  }
#endif
  masked_inverse_s3(mf_p, mf);
  masked_poly_map_mod3_to_modq(&x, mf);


  
  masked_inverse_Sq(&mf_q, &x);

  for(int i=0; i < N_SHARES; ++i) for(int j=0; j < NTRU_N; ++j) mg.shares[i].coeffs[j] = (3*mg.shares[i].coeffs[j])%NTRU_Q;

  sec_Rq_mul(&y, &mg, &mf_q);

  poly_refresh(&y, NTRU_Q);
  
  for(int j = 0; j < NTRU_N; ++j) h.coeffs[j] = y.shares[0].coeffs[j];
  for(int i = 1; i < N_SHARES; ++i) for(int j = 0; j < NTRU_N; ++j) h.coeffs[j] = (h.coeffs[j] + y.shares[i].coeffs[j])%NTRU_Q;
  
  poly_Rq_inv(h_q, &h);
  poly_mod_q_Phi_n(h_q);


  poly_Rq_sum_zero_tobytes(pk, &h);

  uint32_t r;
  for(int i=0; i < N_SHARES*32; i += 4){
    r = rand32();
    *((uint32_t*)(msk_seed + i)) = r;
  }

}




