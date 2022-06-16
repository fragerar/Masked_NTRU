#include "masking.h"
#include "owcpa.h"
#include "poly.h"
#include "sample.h"

static int owcpa_check_ciphertext(const unsigned char *ciphertext)
{
  /* A ciphertext is log2(q)*(n-1) bits packed into bytes.  */
  /* Check that any unused bits of the final byte are zero. */

  uint16_t t = 0;

  t = ciphertext[NTRU_CIPHERTEXTBYTES-1];
  t &= 0xff << (8-(7 & (NTRU_LOGQ*NTRU_PACK_DEG)));

  /* We have 0 <= t < 256 */
  /* Return 0 on success (t=0), 1 on failure */
  return (int) (1&((~t + 1) >> 15));
}

static int owcpa_check_r(const poly *r)
{
  /* A valid r has coefficients in {0,1,q-1} and has r[N-1] = 0 */
  /* Note: We may assume that 0 <= r[i] <= q-1 for all i        */

  int i;
  uint32_t t = 0;
  uint16_t c;
  for(i=0; i<NTRU_N-1; i++)
  {
    c = r->coeffs[i];
    t |= (c + 1) & (NTRU_Q-4);  /* 0 iff c is in {-1,0,1,2} */
    t |= (c + 2) & 4;  /* 1 if c = 2, 0 if c is in {-1,0,1} */
  }
  t |= r->coeffs[NTRU_N-1]; /* Coefficient n-1 must be zero */

  /* We have 0 <= t < 2^16. */
  /* Return 0 on success (t=0), 1 on failure */
  return (int) (1&((~t + 1) >> 31));
}

#ifdef NTRU_HPS
static int owcpa_check_m(const poly *m)
{
  /* Check that m is in message space, i.e.                  */
  /*  (1)  |{i : m[i] = 1}| = |{i : m[i] = 2}|, and          */
  /*  (2)  |{i : m[i] != 0}| = NTRU_WEIGHT.                  */
  /* Note: We may assume that m has coefficients in {0,1,2}. */

  int i;
  uint32_t t = 0;
  uint16_t ps = 0;
  uint16_t ms = 0;
  for(i=0; i<NTRU_N; i++)
  {
    ps += m->coeffs[i] & 1;
    ms += m->coeffs[i] & 2;
  }
  t |= ps ^ (ms >> 1);   /* 0 if (1) holds */
  t |= ms ^ NTRU_WEIGHT; /* 0 if (1) and (2) hold */

  /* We have 0 <= t < 2^16. */
  /* Return 0 on success (t=0), 1 on failure */
  return (int) (1&((~t + 1) >> 31));
}
#endif

void masked_owcpa_enc(unsigned char *c,
               const masked_poly *mr,
               const masked_poly *mm,
               const unsigned char *pk)
{
  int i;
  poly x1, x2;
  poly *h = &x1, *liftm = &x1, *ct = &x2;

  masked_poly x3, x4;
  masked_poly *mct = &x3, *mliftm = &x4;

  poly_Rq_sum_zero_frombytes(h, pk);
  
  halfmasked_poly_Rq_mul(mct, h, mr);
  
  masked_poly_lift(mliftm, mm);
  masked_poly_Rq_add(mct, mct, mliftm);

  unmask_arith_masked_poly(ct, mct, NTRU_Q);
  
  poly_Rq_sum_zero_tobytes(c, ct);
}



int masked_owcpa_dec(unsigned char *rm,
              const unsigned char *ciphertext,
              const masked_poly* mf, const masked_poly* mf_p, const poly* invh)
{
  int i;
  int fail, rm_fail;
  Masked mfail_r;

  uint8_t packed_r[NTRU_PACK_TRINARY_BYTES*(N_SHARES)], packed_m[NTRU_PACK_TRINARY_BYTES*(N_SHARES)];

  poly x1, x2, x3, x4;
  masked_poly mf_q, mcf, mm, mm_lift, mb, mr;

  poly *c = &x1, *f = &x2, *cf = &x3;
  poly *msg_f = &x2, *finv3 = &x3, *m = &x4;
  poly *liftm = &x2, *r = &x4;
  poly *b = &x1;

  poly_Rq_sum_zero_frombytes(c, ciphertext);

  masked_poly_map_mod3_to_modq(&mf_q, mf);



  halfmasked_poly_Rq_mul(&mcf, c, &mf_q); 
  masked_poly_mod3_reduce(&mcf, &mcf);
  sec_S3_mul(&mm, &mcf, mf_p);

  masked_pack_S3(packed_m, &mm);

  fail = 0;

  /* Check that the unused bits of the last byte of the ciphertext are zero */
  fail |= owcpa_check_ciphertext(ciphertext);

  /* For the IND-CCA2 KEM we must ensure that c = Enc(h, (r,m)).             */
  /* We can avoid re-computing r*h + Lift(m) as long as we check that        */
  /* r (defined as b/h mod (q, Phi_n)) and m are in the message space.       */
  /* (m can take any value in S3 in NTRU_HRSS) */
#ifdef NTRU_HPS
  Masked mfail, mfail_message;
  poly_masked_check_message_space(&mfail_message, &mm);
#endif

  /* b = c - Lift(m) mod (q, x^n - 1) */
  masked_poly_lift(&mm_lift, &mm);
  for(i=0; i < NTRU_N; i++){
    mb.shares[0].coeffs[i] = (c->coeffs[i] + NTRU_Q - mm_lift.shares[0].coeffs[i])%NTRU_Q; 
    for(int j=1; j < N_SHARES; ++j)
      mb.shares[j].coeffs[i] = NTRU_Q - mm_lift.shares[j].coeffs[i];
  }
  
  /* r = b / h mod (q, Phi_n) */
  halfmasked_poly_Sq_mul(&mr, invh, &mb);
  pack_s3_from_Rq(packed_r, &mr);

  /* NOTE: Our definition of r as b/h mod (q, Phi_n) follows Figure 4 of     */
  /*   [Sch18] https://eprint.iacr.org/2018/1174/20181203:032458.            */
  /* This differs from Figure 10 of Saito--Xagawa--Yamakawa                  */
  /*   [SXY17] https://eprint.iacr.org/2017/1005/20180516:055500             */
  /* where r gets a final reduction modulo p.                                */
  /* We need this change to use Proposition 1 of [Sch18].                    */

  /* Proposition 1 of [Sch18] shows that re-encryption with (r,m) yields c.  */
  /* if and only if fail==0 after the following call to owcpa_check_r        */
  /* The procedure given in Fig. 8 of [Sch18] can be skipped because we have */
  /* c(1) = 0 due to the use of poly_Rq_sum_zero_{to,from}bytes.             */
  poly_masked_ternary_check(&mfail_r, &mr);



  for(int i=0; i < NTRU_PACK_TRINARY_BYTES; ++i){
    for(int j=0; j < N_SHARES; ++j){
      rm[i                          + j*2*NTRU_PACK_TRINARY_BYTES] = packed_r[i + j*NTRU_PACK_TRINARY_BYTES];
      rm[i+ NTRU_PACK_TRINARY_BYTES + j*2*NTRU_PACK_TRINARY_BYTES] = packed_m[i + j*NTRU_PACK_TRINARY_BYTES];
    }
  }
  #ifdef NTRU_HPS
  sec_and(&mfail, &mfail_message, &mfail_r, 1);
  boolean_refresh(&mfail, 1);
  rm_fail = mfail.shares[0]; 
  for(i=1; i < N_SHARES; ++i) rm_fail ^= mfail.shares[i];
  #endif
  #ifdef NTRU_HRSS
  rm_fail = mfail_r.shares[0]; 
  for(i=1; i < N_SHARES; ++i) rm_fail ^= mfail_r.shares[i];
  #endif
  rm_fail = !rm_fail;
  return fail | (rm_fail);  



}