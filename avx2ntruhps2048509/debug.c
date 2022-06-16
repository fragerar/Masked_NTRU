#include <stdio.h>
#include <stdlib.h>
#include "masking.h"
#include "random.h"


void print_arith_masked(const Masked* x, int q){
  printf(" (");
  int t=0;
  for(int i=0; i < MASKING_ORDER; ++i){
    printf("%u, ",x->shares[i]);
    t += x->shares[i];
  }
  t += x->shares[MASKING_ORDER];
  printf("%u) = %u = %u mod %u\n", x->shares[MASKING_ORDER], t, t%q, q);

}
 

void print_bool_masked(Masked* y){
  int t=0;
  printf(" (");
  for(int i=0; i < MASKING_ORDER; ++i){
    printf("%u, ",y->shares[i]);
    t ^= y->shares[i];
  }
  t ^= y->shares[MASKING_ORDER];
  printf("%i) = (", y->shares[MASKING_ORDER]);
  for(int i=0; i < MASKING_ORDER; ++i){
    printf("0x%X, ", y->shares[i]);
  }
  printf("0x%x) = %u\n", y->shares[MASKING_ORDER], t);
}

void print_unmasked_bool(Masked* y){
  int t=0;
  for(int i=0; i < MASKING_ORDER; ++i){
    t ^= y->shares[i];
  }
  t ^= y->shares[MASKING_ORDER];
  printf(" %u", t); 
}


void print_poly(const poly* p){
  //for(int i=0; i < NTRU_N; ++i) printf("%u ", p->coeffs[i]);
  for(int i=0; i < 20; ++i) printf("%u ", p->coeffs[i]%NTRU_Q);
  printf("\n");
}

void full_print_poly(const poly* p){
  for(int i=0; i < NTRU_N; ++i) printf("%u ", p->coeffs[i]%NTRU_Q);
  
  printf("\n");
}

void full_print_arith_masked_poly(const masked_poly *p, int q){
  uint16_t v;
  Masked t;
  for(int i=0; i < 20; ++i){
    for(int j=0; j < N_SHARES; ++j) t.shares[j] = p->shares[j].coeffs[i];
    print_arith_masked(&t, q);
  }
}

void print_arith_masked_poly(const masked_poly *p, int q){
  uint16_t v;
  for(int i=0; i < 20; ++i){
    v = 0;
    for(int j=0; j < N_SHARES; ++j) v = (v + (p->shares[j]).coeffs[i])%q;
    printf("%u ", v);
  }
  printf("\n");
}

void print_arith_masked_poly_all(const masked_poly *p, int q){
  uint16_t v;
  for(int i=0; i < NTRU_N; ++i){
    v = 0;
    for(int j=0; j < N_SHARES; ++j) v = (v + (p->shares[j]).coeffs[i])%q;
    printf("%u ", v);
  }
  printf("\n");
}

void masked_arith_val(Masked* x, const uint16_t v, int q){
  uint16_t r;

  x->shares[0] = v;
  for(int i=1; i < N_SHARES; ++i){
    r = rand()%q;
    x->shares[i] = r;
    x->shares[0] = (x->shares[0] + q - r)%q;
  }

}

void mask_arith_poly(masked_poly* mp, const poly* p, int q){
  uint16_t v;
  for(int i=0; i < NTRU_N; ++i){ 
    (mp->shares[0]).coeffs[i] = p->coeffs[i];
    for(int j=1; j < N_SHARES; ++j){
      v = rand()%q;
      (mp->shares[0]).coeffs[i] = ((mp->shares[0]).coeffs[i] + q - v)%q; 
      (mp->shares[j]).coeffs[i] = v; 
    } 
  }
}


void unmask_arith_masked_poly(poly *p, const masked_poly* mp, int q){
  uint16_t v;
  for(int i=0; i < NTRU_N; ++i){
    v = 0;
    for(int j=0; j < N_SHARES; ++j) v = (v + (mp->shares[j]).coeffs[i])%q;
    p->coeffs[i] = v;
  }
}


void print_bytes(const uint8_t* bs, int size){
  for(int i=0; i < size; ++i) printf("%X ", bs[i]);
  printf("\n");
}

void mask_bitstring(uint8_t* m_bs, uint8_t* bs, int SIZE){
  uint8_t r;
  for(int i=0; i < SIZE; ++i) m_bs[i] = bs[i];
  for(int i=1; i < N_SHARES; ++i){
    for(int j=0; j < SIZE; ++j) {
      r = (uint8_t) rand16();
      m_bs[i*SIZE + j]  = r;
      m_bs[j]          ^= r;
    }
  }
}

void unmask_bitstring(uint8_t* bs, uint8_t* m_bs, int SIZE){
  for(int i=0; i < SIZE; ++i) bs[i] = m_bs[i];
  for(int i=1; i < N_SHARES; ++i){
    for(int j=0; j < SIZE; ++j) bs[j] ^= m_bs[i*SIZE+j];
  }
}

void print_masked_bs(uint8_t* m_bs, int SIZE){
  uint8_t bs[SIZE*N_SHARES];
  unmask_bitstring(bs, m_bs, SIZE);
  print_bytes(bs, SIZE);
}


