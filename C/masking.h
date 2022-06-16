#ifndef MASKING_H
#define MASKING_H

#include "params.h"
#include "poly.h"

#include <stddef.h>
#include <stdint.h>

#ifndef MASKING_ORDER
#define MASKING_ORDER 1
#endif
#define N_SHARES (MASKING_ORDER+1)


typedef struct Masked {
  uint16_t shares[MASKING_ORDER+1];
} Masked;

typedef struct{
  poly shares[N_SHARES];
} masked_poly;

void masked_owcpa_enc(unsigned char *c,
               const masked_poly *mr,
               const masked_poly *mm,
               const unsigned char *pk);


int masked_owcpa_dec(unsigned char *rm,
              const unsigned char *ciphertext,
              const masked_poly* mf, const masked_poly* mf_p, const poly* h);


int masked_crypto_kem_enc(unsigned char *c, unsigned char *k, const unsigned char *pk);

int masked_crypto_kem_dec(unsigned char *k, const unsigned char *c, const masked_poly* mf, const masked_poly* mf_p, const poly* invh, const uint8_t* msk_seed);


void masked_owcpa_keypair(masked_poly* mf, masked_poly* mf_p, poly* h_q, uint8_t* pk, uint8_t* msk_seed);




void halfmasked_poly_Rq_mul(masked_poly *r, const poly *a, const masked_poly *b);
void halfmasked_poly_Sq_mul(masked_poly *r, const poly *a, const masked_poly *b);
void halfmasked_poly_Rq_add(masked_poly *a, const poly *b);
void masked_poly_Rq_add(masked_poly* c, const masked_poly* a, const masked_poly* b);
void masked_poly_lift(masked_poly* lift_p, const masked_poly* p);
void masked_poly_lift_hrss(masked_poly* lift_p, const masked_poly* p);


void masked_map_mod3_to_modq(Masked* y, const Masked* x, int q); // maps {0,1,2} to {0,1,q-1}
void masked_convert_mod3_to_modq(Masked* y, const Masked* x, int q); // maps {0,1,2} mod 3 to {0,1,2} mod q
void masked_poly_map_mod3_to_modq(masked_poly *b, const masked_poly *a);
void masked_poly_convert_mod3_to_modq(masked_poly *b, const masked_poly *a);

void convert_1bit_B2A_modq(Masked* y, const Masked* x, int q);
void masked_poly_map_mod2_to_modq(masked_poly *b, const masked_poly *a);
void masked_mod3_reduce(Masked* y, const Masked* x);
void masked_poly_mod3_reduce(masked_poly *b, const masked_poly *a);

void masked_pack_S3(uint8_t* masked_bitstring, masked_poly* mp);
void masked_map_modq_2Z3(Masked* y, const Masked* x);
void pack_s3_from_Rq(uint8_t* masked_bitstring, const masked_poly* mp);

void masked_ternary(masked_poly *mp);
void masked_ternary_plus(masked_poly* mp);
void masked_fixed_type(masked_poly* mp, unsigned q);
void uniform_ternary(poly* p);
void uniform_binary(poly* p);
void poly_masked_ternary_check(Masked* b, masked_poly* a);
void poly_masked_check_message_space(Masked* b, const masked_poly* m);

void poly_mod_2_Phi_n(poly *r);
void poly_S2_mul(poly *r, const poly *a, const poly *b);

void sec_Rq_mul(masked_poly* c, const masked_poly* a, const masked_poly* b);
void sec_Sq_mul(masked_poly* c, const masked_poly* a, const masked_poly* b);
void sec_S3_mul(masked_poly* c, const masked_poly* a, const masked_poly* b);
void sec_S2_mul(masked_poly* c, const masked_poly* a, const masked_poly* b);


void linear_arithmetic_refresh(Masked* x, unsigned q);
void poly_refresh(masked_poly* x, unsigned q); 
void sec_and(Masked* res, Masked* x, Masked* y, int k);
void boolean_refresh(Masked* x, unsigned k);


void cube_iter_42(poly* x, poly* y);
void masked_inverse_s3(masked_poly* y, masked_poly* x);
void masked_inverse_s2(masked_poly* y, masked_poly* x);
void masked_inverse_Sq(masked_poly* y, masked_poly* a);


// CGV14
void convert_A2B_CGV14(Masked* y, const Masked* x, unsigned k);
void convert_A2B_CGV14_cadd_1(Masked* y, const Masked* x, unsigned k);

// ------------------- DEBUG ------------------------------
void print_arith_masked_poly(const masked_poly *p, int q);
void print_poly(const poly* p);
void print_unmasked_bool(Masked* y);
void mask_arith_poly(masked_poly* mp, const poly* p, int q);
void unmask_arith_masked_poly(poly *p, const masked_poly* mp, int q);
void print_bytes(const uint8_t* bs, int size);
void print_arith_masked(const Masked* x, int q);
void print_bool_masked(Masked* y);
void full_print_arith_masked_poly(const masked_poly *p, int q);
void unmask_bitstring(uint8_t* bs, uint8_t* m_bs, int SIZE);
void print_masked_bs(uint8_t* m_bs, int SIZE);
void full_print_poly(const poly* p);
void print_arith_masked_poly_all(const masked_poly *p, int q);
void masked_arith_val(Masked* x, const uint16_t v, int q);
void mask_bitstring(uint8_t* m_bs, uint8_t* bs, int SIZE);

#endif