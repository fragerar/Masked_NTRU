#include <stdint.h>

#include "masking.h"
#include "poly.h"
#include "random.h"



/* ================================================== REFRESHING =================================================== */

void linear_arithmetic_refresh(Masked* x, unsigned q){
  int r;
  for(int i=0; i< N_SHARES-1; ++i){

    r = (int) rand32()%q;     
    x->shares[i] = (x->shares[i] + r)%q;
    x->shares[N_SHARES-1] = (x->shares[N_SHARES-1] - r + q)%q;
  }
}

void boolean_refresh(Masked* x, unsigned k){
  int r;
  for(int i=0; i< MASKING_ORDER+1; ++i){
    for(int j=i+1; j < MASKING_ORDER+1; ++j){
      r = (int) rand32() & ((1<<k)-1);
      x->shares[i] = (x->shares[i] ^ r);
      x->shares[j] = (x->shares[j] ^ r);
    }
  }
}


void arithmetic_refresh(Masked* x, unsigned q){
  uint16_t r;

  for(int i=0; i < MASKING_ORDER+1; ++i){
    for(int j=i+1; j < MASKING_ORDER+1; ++j){
      r = rand16();     
      x->shares[i] = (x->shares[i] + r)%q;
      x->shares[j] = (x->shares[j] - r + q)%q;
    }
  }
}

void poly_arithmetic_refresh(masked_poly* mp, unsigned q){
  uint16_t r;

  for(int k=0; k < NTRU_N; ++k){
    for(int i=0; i < MASKING_ORDER+1; ++i){
      for(int j=i+1; j < MASKING_ORDER+1; ++j){
        r = rand16()%q;     
        mp->shares[i].coeffs[k] = (mp->shares[i].coeffs[k] + r)%q;
        mp->shares[j].coeffs[k] = (mp->shares[j].coeffs[k] - r + q)%q;
      }
    }
  }
}

void poly_refresh(masked_poly* x, unsigned q){

  Masked y;
  for(int i=0; i < NTRU_N; ++i){
    for(int j=0; j < N_SHARES; ++j) y.shares[j] = x->shares[j].coeffs[i];
    linear_arithmetic_refresh(&y, q);
    for(int j=0; j < N_SHARES; ++j) x->shares[j].coeffs[i] = y.shares[j];
  }

}
/* ================================================================================================================= */







/* ==================================================== BASICS ===================================================== */


void sec_and(Masked* res, Masked* x, Masked* y, int k){

#if MASKING_ORDER == 1
    uint16_t u = rand16()&((1<<k)-1);
    uint16_t z;
    z = u ^ (x->shares[0] & y->shares[0]);
    z = z ^ (x->shares[0] & y->shares[1]);
    z = z ^ (x->shares[1] & y->shares[0]);
    z = z ^ (x->shares[1] & y->shares[1]);
    res->shares[0] = z;
    res->shares[1] = u;

#else
    Masked r;
    uint16_t i, j, z_ij, z_ji;
    for(i=0; i < MASKING_ORDER + 1; ++i) r.shares[i] = x->shares[i] & y->shares[i];
    for(i=0; i < MASKING_ORDER + 1; ++i)
        for(j=i+1; j < MASKING_ORDER + 1; ++j){
            z_ij  = rand16()&((1<<k)-1);
            z_ji  = (x->shares[i] & y->shares[j]) ^ z_ij;
            z_ji ^= (x->shares[j] & y->shares[i]);
            r.shares[i] ^= z_ij;
            r.shares[j] ^= z_ji;            
        }
    for(i=0; i < MASKING_ORDER + 1; ++i) res->shares[i] = r.shares[i];
#endif

}


void sec_mult(Masked* c, const Masked* a, const Masked* b, unsigned q){
  uint16_t r, t;
 
  for(int i=0; i < MASKING_ORDER+1; ++i) c->shares[i] = (a->shares[i]*b->shares[i])%q;

  for(int i=0; i < MASKING_ORDER+1; ++i){

    for(int j=i+1; j < MASKING_ORDER+1; ++j){
      r = rand16();
      t = ((r+a->shares[i]*b->shares[j])+a->shares[j]*b->shares[i])%q;
      c->shares[i] = (c->shares[i] + q - r)%q; 
      c->shares[j] = (c->shares[j] + t)%q; 
    }
  }
}

/* ================================================================================================================= */







/* ================================================ POLY ARITHMETIC ================================================ */


void halfmasked_poly_Rq_mul(masked_poly *r, const poly *a, const masked_poly *b){
  for(int i=0; i < N_SHARES; ++i){
    poly_Rq_mul(&(r->shares[i]), a, &(b->shares[i]));
  }
}

void halfmasked_poly_Sq_mul(masked_poly *r, const poly *a, const masked_poly *b){
  for(int i=0; i < N_SHARES; ++i){
    poly_Sq_mul(&(r->shares[i]), a, &(b->shares[i]));
  }
}

void halfmasked_poly_Rq_add(masked_poly *a, const poly *b){
  for(int i=0; i < NTRU_N; ++i) 
    (a->shares[0]).coeffs[i] = ((a->shares[0]).coeffs[i] + b->coeffs[i])%NTRU_Q; 
}

void masked_poly_Rq_add(masked_poly* c, const masked_poly* a, const masked_poly* b){
  for(int i=0; i < N_SHARES; ++i){
    for(int j=0; j < NTRU_N; ++j) c->shares[i].coeffs[j] = (a->shares[i].coeffs[j] + b->shares[i].coeffs[j])%NTRU_Q; 
  }
}



#ifdef NTRU_HRSS
static void poly_divide_by_x_minus_one(poly *r, const poly *a)
{
  int i;
  poly b;
  uint16_t t, zj;
  t = 3 - (NTRU_N % 3);
  b.coeffs[0] = a->coeffs[0] * (2-t) + a->coeffs[1] * 0 + a->coeffs[2] * t;
  b.coeffs[1] = a->coeffs[1] * (2-t) + a->coeffs[2] * 0;
  b.coeffs[2] = a->coeffs[2] * (2-t);

  zj = 0; /* z[1] */
  for(i=3; i<NTRU_N; i++)
  {
    b.coeffs[0] += a->coeffs[i] * (zj + 2*t);
    b.coeffs[1] += a->coeffs[i] * (zj + t);
    b.coeffs[2] += a->coeffs[i] * zj;
    zj = (zj + t) % 3;
  }
  b.coeffs[1] += a->coeffs[0] * (zj + t);
  b.coeffs[2] += a->coeffs[0] * zj;
  b.coeffs[2] += a->coeffs[1] * (zj + t);

  b.coeffs[0] = b.coeffs[0];
  b.coeffs[1] = b.coeffs[1];
  b.coeffs[2] = b.coeffs[2];
  for(i=3; i<NTRU_N; i++) {
    b.coeffs[i] = b.coeffs[i-3] + 2*(a->coeffs[i] + a->coeffs[i-1] + a->coeffs[i-2]);
  }

  poly_mod_3_Phi_n(&b);

  for(int i=0; i < NTRU_N; ++i) r->coeffs[i] = b.coeffs[i];
}
#endif

void poly_rand(poly* x){
  for(int i=0; i < NTRU_N; ++i) x->coeffs[i] = rand16()%NTRU_Q;
}

void poly_add_modq(poly* c, const poly* a, const poly* b, unsigned q){
  for(int i=0; i < NTRU_N; ++i) c->coeffs[i] = (a->coeffs[i] + b->coeffs[i])%q; 
}

void poly_sub_modq(poly* c, const poly* a, const poly* b, unsigned q){
  for(int i=0; i < NTRU_N; ++i) c->coeffs[i] = (a->coeffs[i] + q - b->coeffs[i])%q; 
}


void sec_Rq_mul(masked_poly* c, const masked_poly* a, const masked_poly* b){
  poly r, t, x;
 
  for(int i=0; i < MASKING_ORDER+1; ++i) poly_Rq_mul(&(c->shares[i]), &(a->shares[i]), &(b->shares[i]));

  for(int i=0; i < MASKING_ORDER+1; ++i){

    for(int j=i+1; j < MASKING_ORDER+1; ++j){
      poly_rand(&r);
      poly_Rq_mul(&x, &(a->shares[i]), &(b->shares[j]));
      poly_add_modq(&t, &r, &x, NTRU_Q);
      poly_Rq_mul(&x, &(a->shares[j]), &(b->shares[i]));
      poly_add_modq(&t, &t, &x, NTRU_Q);
      poly_sub_modq(&(c->shares[i]), &(c->shares[i]), &r, NTRU_Q);
      poly_add_modq(&(c->shares[j]), &(c->shares[j]), &t, NTRU_Q);
    }
  }
} 

void sec_Sq_mul(masked_poly* c, const masked_poly* a, const masked_poly* b){
  poly r, t, x;
 
  for(int i=0; i < MASKING_ORDER+1; ++i) poly_Sq_mul(&(c->shares[i]), &(a->shares[i]), &(b->shares[i]));

  for(int i=0; i < MASKING_ORDER+1; ++i){

    for(int j=i+1; j < MASKING_ORDER+1; ++j){
      poly_rand(&r);
      poly_Sq_mul(&x, &(a->shares[i]), &(b->shares[j]));
      poly_add_modq(&t, &r, &x, NTRU_Q);
      poly_Sq_mul(&x, &(a->shares[j]), &(b->shares[i]));
      poly_add_modq(&t, &t, &x, NTRU_Q);
      poly_sub_modq(&(c->shares[i]), &(c->shares[i]), &r, NTRU_Q);
      poly_add_modq(&(c->shares[j]), &(c->shares[j]), &t, NTRU_Q);
    }
  }
} 


void poly_mod_2_Phi_n(poly *r)
{
  int i;
  for(i=0; i <NTRU_N; i++)
    r->coeffs[i] = (r->coeffs[i] ^ r->coeffs[NTRU_N-1])&1;
}

void poly_S2_mul(poly *r, const poly *a, const poly *b)
{
  poly_Rq_mul(r, a, b);
  poly_mod_2_Phi_n(r);
}


void sec_S2_mul(masked_poly* c, const masked_poly* a, const masked_poly* b){
  poly r, t, x;
 
  for(int i=0; i < MASKING_ORDER+1; ++i) poly_S2_mul(&(c->shares[i]), &(a->shares[i]), &(b->shares[i]));

  for(int i=0; i < MASKING_ORDER+1; ++i){

    for(int j=i+1; j < MASKING_ORDER+1; ++j){
      uniform_binary(&r);
      poly_S2_mul(&x, &(a->shares[i]), &(b->shares[j]));
      poly_add_modq(&t, &r, &x, 2);
      poly_S2_mul(&x, &(a->shares[j]), &(b->shares[i]));
      poly_add_modq(&t, &t, &x, 2);
      poly_sub_modq(&(c->shares[i]), &(c->shares[i]), &r, 2);
      poly_add_modq(&(c->shares[j]), &(c->shares[j]), &t, 2);
    }
  }
} 




void sec_S3_mul(masked_poly* c, const masked_poly* a, const masked_poly* b){
  poly r, t, x;
 
  for(int i=0; i < MASKING_ORDER+1; ++i) poly_S3_mul(&(c->shares[i]), &(a->shares[i]), &(b->shares[i]));

  for(int i=0; i < MASKING_ORDER+1; ++i){

    for(int j=i+1; j < MASKING_ORDER+1; ++j){
      uniform_ternary(&r);
      poly_S3_mul(&x, &(a->shares[i]), &(b->shares[j]));
      poly_add_modq(&t, &r, &x, 3);
      poly_S3_mul(&x, &(a->shares[j]), &(b->shares[i]));
      poly_add_modq(&t, &t, &x, 3);
      poly_sub_modq(&(c->shares[i]), &(c->shares[i]), &r, 3);
      poly_add_modq(&(c->shares[j]), &(c->shares[j]), &t, 3);
    }
  }
} 

/* ================================================================================================================= */







/* ====================================================== LIFT ===================================================== */

#ifdef NTRU_HPS
void masked_poly_lift(masked_poly* lift_p, const masked_poly* p){
  masked_poly_map_mod3_to_modq(lift_p, p);
}
#endif

#ifdef NTRU_HRSS
void masked_poly_lift(masked_poly* lift_p, const masked_poly* p){
  masked_poly b;
  for(int i=0; i < N_SHARES; ++i) poly_divide_by_x_minus_one(&(lift_p->shares[i]), &(p->shares[i]));
  masked_poly_map_mod3_to_modq(&b, lift_p);

  for(int k=0; k < N_SHARES; ++k)
  {
    lift_p->shares[k].coeffs[0] = -(b.shares[k].coeffs[0]);
    for(int i=0; i<NTRU_N-1; i++) {
      lift_p->shares[k].coeffs[i+1] = b.shares[k].coeffs[i] - b.shares[k].coeffs[i+1];
    }
  }
}
#endif

/* ================================================================================================================= */






/* ================================================== CONVERSIONS ================================================== */

void masked_map_mod3_to_modq(Masked* y, const Masked* x, int q){
  // {0, 1, 2} -> {0, 1, q-1}
  Masked T[3];
  Masked T_p[3];

  T[0].shares[0] = 0;
  for(int i=1; i < N_SHARES; ++i) T[0].shares[i] = 0;

  T[1].shares[0] = 1;
  for(int i=1; i < N_SHARES; ++i) T[1].shares[i] = 0;

  T[2].shares[0] = q-1;
  for(int i=1; i < N_SHARES; ++i) T[2].shares[i] = 0;

  for(int i=0; i < N_SHARES-1; ++i){
    for(int u=0; u < 3; ++u){
      for(int j=0; j < N_SHARES; ++j){
        T_p[u].shares[j] = T[(u+(x->shares[i]))%3].shares[j]; 
      }
    }
    for(int u=0; u < 3; ++u){
      linear_arithmetic_refresh(&(T_p[u]), q); 
      for(int i=0; i < N_SHARES; ++i) T[u].shares[i] = T_p[u].shares[i];
    }
  }
  for(int i=0; i < N_SHARES; ++i) y->shares[i] = T[x->shares[MASKING_ORDER]].shares[i]; 
  linear_arithmetic_refresh(y, q);

}

void masked_poly_map_mod3_to_modq(masked_poly *b, const masked_poly *a){
  // Maps Z_3[x]/f to Z_q[x]/f using the mapping {0, 1, 2} -> {0, 1, q-1}
  Masked x,y;
  for(int i=0; i < NTRU_N; ++i){
    for(int j=0; j < N_SHARES; ++j) x.shares[j] = a->shares[j].coeffs[i];
    masked_map_mod3_to_modq(&y, &x, NTRU_Q);
    for(int j=0; j < N_SHARES; ++j) b->shares[j].coeffs[i] = y.shares[j];
  }
}

void convert_1bit_B2A_modq(Masked* y, const Masked* x, int q){
  Masked T[2];
  Masked T_p[2];
  
  T[0].shares[0] = 0;
  for(int i=1; i < N_SHARES; ++i) T[0].shares[i] = 0;

  T[1].shares[0] = 1;
  for(int i=1; i < N_SHARES; ++i) T[1].shares[i] = 0;

  
  for(int i=0; i < N_SHARES-1; ++i){
    for(int u=0; u < 2; ++u){
      for(int j=0; j < N_SHARES; ++j){
        T_p[u].shares[j] = T[(u^(x->shares[i]))].shares[j]; 
      }
    }
    for(int u=0; u < 2; ++u){
      linear_arithmetic_refresh(&(T_p[u]), q); 
      for(int i=0; i < N_SHARES; ++i) T[u].shares[i] = T_p[u].shares[i];
    }
  }
  for(int i=0; i < N_SHARES; ++i) y->shares[i] = T[x->shares[MASKING_ORDER]].shares[i]; 
  linear_arithmetic_refresh(y, q); 
}



void masked_poly_map_mod2_to_modq(masked_poly *b, const masked_poly *a){
  // Maps Z_2[x]/f to Z_q[x]/f using the mapping {0, 1} mod 2 -> {0, 1} mod q
  Masked x,y;
  for(int i=0; i < NTRU_N; ++i){
    for(int j=0; j < N_SHARES; ++j) x.shares[j] = a->shares[j].coeffs[i];
    convert_1bit_B2A_modq(&y, &x, NTRU_Q);
    for(int j=0; j < N_SHARES; ++j) b->shares[j].coeffs[i] = y.shares[j];
  }
}



void masked_convert_mod3_to_modq(Masked* y, const Masked* x, int q){
  // {0, 1, 2} mod 3 -> {0, 1, 2} mod q
  Masked T[3];
  Masked T_p[3];

  T[0].shares[0] = 0;
  for(int i=1; i < N_SHARES; ++i) T[0].shares[i] = 0;

  T[1].shares[0] = 1;
  for(int i=1; i < N_SHARES; ++i) T[1].shares[i] = 0;

  T[2].shares[0] = 2;
  for(int i=1; i < N_SHARES; ++i) T[2].shares[i] = 0;

  for(int i=0; i < N_SHARES-1; ++i){
    for(int u=0; u < 3; ++u){
      for(int j=0; j < N_SHARES; ++j){
        T_p[u].shares[j] = T[(u+(x->shares[i]))%3].shares[j]; 
      }
    }
    for(int u=0; u < 3; ++u){
      linear_arithmetic_refresh(&(T_p[u]), q); 
      for(int i=0; i < N_SHARES; ++i) T[u].shares[i] = T_p[u].shares[i];
    }
  }
  for(int i=0; i < N_SHARES; ++i) y->shares[i] = T[x->shares[MASKING_ORDER]].shares[i]; 
  linear_arithmetic_refresh(y, q);

}


void masked_poly_convert_mod3_to_modq(masked_poly *b, const masked_poly *a){
  Masked x,y;
  for(int i=0; i < NTRU_N; ++i){
    for(int j=0; j < N_SHARES; ++j) x.shares[j] = a->shares[j].coeffs[i];
    masked_convert_mod3_to_modq(&y, &x, NTRU_Q);
    for(int j=0; j < N_SHARES; ++j) b->shares[j].coeffs[i] = y.shares[j];
  }
}/* ================================================================================================================= */







/* ==================================================== PACKING ==================================================== */

void masked_pack_S3(uint8_t* masked_bitstring, masked_poly* mp){
  Masked v[5];
  Masked t, packed, packed_bool;
  int k;
  for(int i=0; i < NTRU_PACK_DEG; i+=5){
    for(int j=0; j < N_SHARES; ++j){
      k = 0;
      while (((i+k) < NTRU_PACK_DEG) && (k < 5)){
        v[k].shares[j] = mp->shares[j].coeffs[i+k];
        ++k;
      }
      while (k < 5){
        v[k].shares[j] = 0;
        ++k;
      }
    }
    masked_convert_mod3_to_modq(&packed, &v[0], 256);
    masked_convert_mod3_to_modq(&t, &v[1], 256);
    for(int j=0; j < N_SHARES; ++j) packed.shares[j] = (packed.shares[j] + 3*t.shares[j])%256;
    masked_convert_mod3_to_modq(&t, &v[2], 256);
    for(int j=0; j < N_SHARES; ++j) packed.shares[j] = (packed.shares[j] + 9*t.shares[j])%256;
    masked_convert_mod3_to_modq(&t, &v[3], 256);
    for(int j=0; j < N_SHARES; ++j) packed.shares[j] = (packed.shares[j] + 27*t.shares[j])%256;
    masked_convert_mod3_to_modq(&t, &v[4], 256);
    for(int j=0; j < N_SHARES; ++j) packed.shares[j] = (packed.shares[j] + 81*t.shares[j])%256;


    convert_A2B_CGV14(&packed_bool, &packed, 8);
    for(int j=0; j < N_SHARES; ++j) masked_bitstring[j*NTRU_PACK_TRINARY_BYTES+i/5] = packed_bool.shares[j];
  }
}


void masked_map_modq_2Z3(Masked* y, const Masked* x){
  //x*(511 + 3*x)
  Masked t, t2;
  for(int i=0; i < N_SHARES; ++i) t.shares[i] = (x->shares[i]&511);
  for(int i=0; i < N_SHARES; ++i) t2.shares[i] = (3*t.shares[i])&511;
  t2.shares[0] = (t2.shares[0] + 511)&511;
  arithmetic_refresh(&t2, 512);
  sec_mult(y, &t, &t2, 512);  
}

void pack_s3_from_Rq(uint8_t* masked_bitstring, const masked_poly* mp){

  Masked v[5], t, res, packed_bool;
  int k;
  for(int i=0; i < NTRU_PACK_DEG; i+=5){
    for(int j=0; j < N_SHARES; ++j){
      k = 0;
      while (((i+k) < NTRU_PACK_DEG) && (k < 5)){
        v[k].shares[j] = mp->shares[j].coeffs[i+k];
        ++k;
      }
      while (k < 5){
        v[k].shares[j] = 0;
        ++k;
      }
    }

    masked_map_modq_2Z3(&t, &v[0]);
    for(int j=0; j < N_SHARES; ++j) res.shares[j] = t.shares[j];

    masked_map_modq_2Z3(&t, &v[1]);
    for(int j=0; j < N_SHARES; ++j) res.shares[j] = (res.shares[j] + 3*t.shares[j])&511;

    masked_map_modq_2Z3(&t, &v[2]);
    for(int j=0; j < N_SHARES; ++j) res.shares[j] = (res.shares[j] + 9*t.shares[j])&511;

    masked_map_modq_2Z3(&t, &v[3]);
    for(int j=0; j < N_SHARES; ++j) res.shares[j] = (res.shares[j] + 27*t.shares[j])&511;

    masked_map_modq_2Z3(&t, &v[4]);
    for(int j=0; j < N_SHARES; ++j) res.shares[j] = (res.shares[j] + 81*t.shares[j])&511;
    
    convert_A2B_CGV14(&packed_bool, &res, 9);
    for(int j=0; j < N_SHARES; ++j) masked_bitstring[j*NTRU_PACK_TRINARY_BYTES+i/5] = (uint8_t)(packed_bool.shares[j]>>1);
    
  }
}
/* ================================================================================================================= */







/* ================================================ REDUCTION MOD 3 ================================================ */


void masked_mod3_reduce(Masked* y, const Masked* x){
  // Masked reduction modulo 3 of the representative in [-q/2, q/2) of an element of Z_q 
  Masked t1, t2, x_bool;
  convert_A2B_CGV14(&x_bool, x, NTRU_LOGQ);
  for(int j=0; j < N_SHARES; ++j) t1.shares[j] =  (x_bool.shares[j])&1;
  convert_1bit_B2A_modq(&t2, &t1, 3);
  for(int j=0; j < N_SHARES; ++j) y->shares[j] = (t2.shares[j])%3;

  for(int i=1; i < NTRU_LOGQ-1; ++i){
    for(int j=0; j < N_SHARES; ++j) t1.shares[j] =  (x_bool.shares[j] >> i)&1;
    convert_1bit_B2A_modq(&t2, &t1, 3);
    for(int j=0; j < N_SHARES; ++j) y->shares[j] = (y->shares[j] + (1<<i)*t2.shares[j])%3;
  }

  for(int j=0; j < N_SHARES; ++j) t1.shares[j] =  (x_bool.shares[j] >> (NTRU_LOGQ-1))&1;
  convert_1bit_B2A_modq(&t2, &t1, 3);
  for(int j=0; j < N_SHARES; ++j) y->shares[j] = (y->shares[j] + (3-((1<<(NTRU_LOGQ-1))%3))*t2.shares[j])%3; 
}


void masked_poly_mod3_reduce(masked_poly *b, const masked_poly *a){
  Masked x, y;
  for(int i=0; i < NTRU_N; ++i){
    for(int j=0; j < N_SHARES; ++j) x.shares[j] = (a->shares[j].coeffs[i])%NTRU_Q;
    masked_mod3_reduce(&y, &x);
    for(int j=0; j < N_SHARES; ++j) b->shares[j].coeffs[i] = y.shares[j];
  }
}

/* ================================================================================================================= */







/* ================================================= ZT AND CHECKS ================================================= */
void boolean_zero_test(Masked* y, const Masked* x, int k, int logk){
  Masked z;
  y->shares[0] = (~x->shares[0]) | ((1<<(1<<logk))-(1<<k));
  for(int i=1; i < MASKING_ORDER+1; ++i) y->shares[i] = x->shares[i];

  for(int i=0; i < logk; ++i){
    for(int j=0; j < MASKING_ORDER+1; ++j) z.shares[j] = y->shares[j] >> (1 << i);
    boolean_refresh(&z, k);
    sec_and(y, &z, y, k); 
  }
  for(int i=0; i < MASKING_ORDER+1; ++i) y->shares[i] = (y->shares[i])&1;

  boolean_refresh(y, k);

}

void masked_ternary_check(Masked* y, Masked* x){
  Masked x_bool, b0, b1, b2;
  convert_A2B_CGV14(&x_bool, x, NTRU_LOGQ);
  
  boolean_zero_test(&b0, &x_bool, NTRU_LOGQ, 4); // Check 0
  x_bool.shares[0] ^= 1;
  boolean_zero_test(&b1, &x_bool, NTRU_LOGQ, 4); // Check 1
  x_bool.shares[0] ^= 1^(NTRU_Q-1);
  boolean_zero_test(&b2, &x_bool, NTRU_LOGQ, 4); // Check -1

  for(int i=0; i < N_SHARES; ++i) y->shares[i] = b0.shares[i] ^ b1.shares[i] ^ b2.shares[i]; 

}

void poly_masked_ternary_check(Masked* b, masked_poly* a){
  Masked x, y;

  b->shares[0] = 1;
  for(int i=1; i < N_SHARES; ++i) b->shares[i] = 0;

  for(int i=0; i < NTRU_N; ++i){
    for(int j=0; j < N_SHARES; ++j) x.shares[j] = a->shares[j].coeffs[i];
    masked_ternary_check(&y, &x);
    sec_and(&x, b, &y, 1);
    for(int j=0; j < N_SHARES; ++j) b->shares[j] = x.shares[j];
  }
}

void sec_square(Masked* y, const Masked* x, int q){
  Masked fresh_x;
  for(int i=0; i < N_SHARES; ++i) fresh_x.shares[i] = x->shares[i];
  arithmetic_refresh(&fresh_x, q);
  sec_mult(y, &fresh_x, x, q);
} 

void masked_arith_equality_test(Masked* b, Masked* x, uint16_t val, unsigned k, unsigned logk){
  // testing if x is a masking of val
  Masked x_bool;
  convert_A2B_CGV14(&x_bool, x, k);
  x_bool.shares[0] ^= val;
  boolean_zero_test(b, &x_bool, k, logk);
}

void poly_masked_check_message_space(Masked* b, const masked_poly* m){
  Masked sum_square, sum, t, x;
  Masked b0, b1;
  masked_poly lift_m;

  masked_poly_lift(&lift_m, m);
  

  for(int i=0; i < N_SHARES; ++i){
    sum_square.shares[i] = 0;
    sum.shares[i]        = 0;
  }

  for(int i=0; i < NTRU_N; ++i){
    for(int j=0; j < N_SHARES; ++j) x.shares[j] = lift_m.shares[j].coeffs[i];
    sec_square(&t, &x, NTRU_Q);
    for(int j=0; j < N_SHARES; ++j){
      sum.shares[j]        = (sum.shares[j]        + x.shares[j])%NTRU_Q;
      sum_square.shares[j] = (sum_square.shares[j] + t.shares[j])%NTRU_Q;
    }
  }

  masked_arith_equality_test(&b0, &sum, 0, NTRU_LOGQ, 4);  
  masked_arith_equality_test(&b1, &sum_square, NTRU_WEIGHT, NTRU_LOGQ, 4);  
  sec_and(b, &b0, &b1, 1);
}

/* ================================================================================================================= */







/* =================================================== SAMPLING ==================================================== */


void uniform_ternary(poly* p){
  uint32_t r;
  uint8_t x;
  int stop=0;

  int i=0;
  while (!stop){
    r = rand32();
    for(int j=0; (j < 16) & (!stop); ++j){
      x = r&0x03;
      r >>= 2;
      if (x != 3) {
        p->coeffs[i] = x;
        i++;
        if (i==NTRU_N-1) stop = 1;
      } 
    } 
  }
  p->coeffs[NTRU_N-1] = 0;
}

void uniform_binary(poly* p){
  // Should mprove randomness usage !
  uint32_t r;
  for(int i=0; i < NTRU_N-1; ++i) p->coeffs[i] = rand16()&1;
  p->coeffs[NTRU_N-1] = 0;
}

void masked_ternary(masked_poly *mp){
  for(int i=0; i < N_SHARES; ++i) uniform_ternary(&(mp->shares[i]));
}

uint16_t random_uniform_number16bits(uint16_t n){
  /* Uniform over [0, n-1] */
  int bound = ((1<<16)/n)*n;
  uint16_t r = rand16();

  while (r >= bound) r = rand16();

  return r%n;

}

#ifdef NTRU_HPS

void masked_fixed_type(masked_poly* mp, unsigned q){
  int x;
  uint16_t temp;
  for (int i =             0; i < NTRU_WEIGHT/2; i++) mp->shares[0].coeffs[i] =  1;
  for (int i = NTRU_WEIGHT/2; i <   NTRU_WEIGHT; i++) mp->shares[0].coeffs[i] =  2;
  for (int i =   NTRU_WEIGHT; i <        NTRU_N; i++) mp->shares[0].coeffs[i] =  0;

  for(int i=0; i < NTRU_N; ++i){
    for(int j=1; j < N_SHARES; ++j) mp->shares[j].coeffs[i] = 0;
  }
  
  for(int k=0; k < N_SHARES; ++k){
    for(int i=NTRU_N-1; i > 0; --i){ // Fisher-Yates shuffle
      x = random_uniform_number16bits(i);
      for(int j=0; j < N_SHARES; ++j){
        temp = mp->shares[j].coeffs[x];
        mp->shares[j].coeffs[x] = mp->shares[j].coeffs[i-1];
        mp->shares[j].coeffs[i-1] = temp;
      }
    }
    poly_arithmetic_refresh(mp, q);
  }
}

#endif


void masked_ternary_plus(masked_poly* mp){
  Masked t, x1, x2, c;
  Masked *coeffA, *coeffB, *swap;
  masked_ternary(mp);
  int s = 0;

  coeffA = &x1; coeffB = &x2;
  
  for(int i=0; i < N_SHARES; ++i) t.shares[i] = 0;

  for(int i=0; i < N_SHARES; ++i) coeffA->shares[i] = mp->shares[i].coeffs[0];
  masked_map_mod3_to_modq(coeffA, coeffA, NTRU_Q);  

  for(int i=1; i < NTRU_N; ++i){
    for(int j=0; j < N_SHARES; ++j) coeffB->shares[j] = mp->shares[j].coeffs[i];
    masked_map_mod3_to_modq(coeffB, coeffB, NTRU_Q);
    sec_mult(&c, coeffA, coeffB, NTRU_Q);
    for(int j=0; j < N_SHARES; ++j) t.shares[j] = (t.shares[j] + c.shares[j])%NTRU_Q;

    swap = coeffA;
    coeffA = coeffB;
    coeffB = swap;

  }
  
  convert_A2B_CGV14(&t, &t, NTRU_LOGQ);

  for(int i=0; i < N_SHARES; ++i) t.shares[i] >>= (NTRU_LOGQ - 1); 
  boolean_refresh(&t, 1);
  for(int i=0; i < N_SHARES; ++i) s ^= t.shares[i];
  s += 1; //If s=0: s=1 else if s=1; s=2 = (-1) mod 3
  for(int i=0; i<NTRU_N; i+=2){
    for(int j=0; j < N_SHARES; ++j)
      mp->shares[j].coeffs[i] = (s*mp->shares[j].coeffs[i])%3; 
  }

}

void masked_sample_rm(masked_poly* mr, masked_poly* mm){
#ifdef NTRU_HPS
  masked_ternary(mr);
  masked_fixed_type(mm, 3);
#endif

#ifdef NTRU_HRSS
  masked_ternary_plus(mr);
  masked_ternary_plus(mm);
#endif
}


/* ================================================================================================================= */










