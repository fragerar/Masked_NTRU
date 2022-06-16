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



static void dummy_unmasked_keygen(masked_poly* mf, masked_poly* mg){

  printf("======== Dummy masked ========\n");

  int i;

  poly test;
  poly x1, x2, x3, x4, x5;

  poly *f=&x1, *g=&x2, *invf_mod3=&x3;
  poly *gf=&x3, *invgf=&x4, *tmp=&x5;
  poly *invh=&x3, *h=&x3;


  unmask_arith_masked_poly(f, mf, 3);
#ifdef NTRU_HRSS
  unmask_arith_masked_poly(g, mg, 3);
#endif

#ifdef NTRU_HPS
  unmask_arith_masked_poly(g, mg, 3);
#endif



  poly_S3_inv(invf_mod3, f);
  printf("mf_p: "); print_poly(invf_mod3);

  // Lift coeffs of f and g from Z_p to Z_q 
  poly_Z3_to_Zq(f);
  poly_Z3_to_Zq(g);
  printf("f: "); print_poly(f);
#ifdef NTRU_HRSS
  // g = 3*(x-1)*g 
  for(i=NTRU_N-1; i>0; i--)
    g->coeffs[i] = 3*(g->coeffs[i-1] - g->coeffs[i]);
  g->coeffs[0] = -(3*g->coeffs[0]);
#endif
  



#ifdef NTRU_HPS
  // g = 3*g 
  for(i=0; i<NTRU_N; i++)
    g->coeffs[i] = 3 * g->coeffs[i];
#endif

  printf("g: "); print_poly(g);

  poly_Rq_mul(gf, g, f);
  poly_Rq_inv(invgf, gf);

  poly_Rq_mul(tmp, invgf, f);
  poly_Sq_mul(invh, tmp, f);

  poly_Rq_mul(tmp, invgf, g);
  poly_Rq_mul(h, tmp, g);

  printf("h: "); print_poly(h);
  printf("=============================\n");


} 


void masked_owcpa_keypair(masked_poly* mf, masked_poly* mf_p, poly* h_q, uint8_t* pk, uint8_t* msk_seed){
  poly h;
  masked_poly mg, mf_q, x, y;

#ifdef NTRU_HPS
  masked_ternary(mf);
  masked_fixed_type(&y, 3);
  //dummy_unmasked_keygen(mf, &y);

  masked_poly_map_mod3_to_modq(&mg, &y);
#endif

#ifdef NTRU_HRSS
  masked_ternary_plus(mf);
  masked_ternary_plus(&mg);

  //dummy_unmasked_keygen(mf, &mg);

  masked_poly_map_mod3_to_modq(&y, &mg);


  for(int j=0; j < N_SHARES; ++j){
    for(int i=NTRU_N-1; i>0; i--)
      mg.shares[j].coeffs[i] = (y.shares[j].coeffs[i-1] + NTRU_Q - y.shares[j].coeffs[i])%NTRU_Q;
    mg.shares[j].coeffs[0] = NTRU_Q-(y.shares[j].coeffs[0]);
  }
#endif
  masked_inverse_s3(mf_p, mf);
  //printf("masked_f_p: "); print_arith_masked_poly(mf_p, 3);
  masked_poly_map_mod3_to_modq(&x, mf);


  //printf("masked_f: "); print_arith_masked_poly(&x, NTRU_Q);
  
  masked_inverse_Sq(&mf_q, &x);

  for(int i=0; i < N_SHARES; ++i) for(int j=0; j < NTRU_N; ++j) mg.shares[i].coeffs[j] = (3*mg.shares[i].coeffs[j])%NTRU_Q;

  //printf("masked_g: "); print_arith_masked_poly(&mg, NTRU_Q);
  sec_Rq_mul(&y, &mg, &mf_q);

  poly_refresh(&y, NTRU_Q);
  
  for(int j = 0; j < NTRU_N; ++j) h.coeffs[j] = y.shares[0].coeffs[j];
  for(int i = 1; i < N_SHARES; ++i) for(int j = 0; j < NTRU_N; ++j) h.coeffs[j] = (h.coeffs[j] + y.shares[i].coeffs[j])%NTRU_Q;
  
  poly_Rq_inv(h_q, &h);
  poly_mod_q_Phi_n(h_q);


  poly_Rq_sum_zero_tobytes(pk, &h);
  //printf("mf: "); print_arith_masked_poly(mf, 3);
  //printf("masked_h: ");print_poly(&h);

  uint32_t r;
  for(int i=0; i < N_SHARES*32; i += 4){
    r = rand32();
    *((uint32_t*)(msk_seed + i)) = r;
  }

}




