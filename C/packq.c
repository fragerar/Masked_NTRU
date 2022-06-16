#include "poly.h"

#if defined(PARAMS_1) || defined(PARAMS_2)
void poly_Sq_tobytes(unsigned char *r, const poly *a)
{
  int i,j;
  uint16_t t[8];

  for(i=0;i<NTRU_PACK_DEG/8;i++)
  {
    for(j=0;j<8;j++)
      t[j] = MODQ(a->coeffs[8*i+j]);

    r[11 * i + 0] = (unsigned char) ( t[0]        & 0xff);
    r[11 * i + 1] = (unsigned char) ((t[0] >>  8) | ((t[1] & 0x1f) << 3));
    r[11 * i + 2] = (unsigned char) ((t[1] >>  5) | ((t[2] & 0x03) << 6));
    r[11 * i + 3] = (unsigned char) ((t[2] >>  2) & 0xff);
    r[11 * i + 4] = (unsigned char) ((t[2] >> 10) | ((t[3] & 0x7f) << 1));
    r[11 * i + 5] = (unsigned char) ((t[3] >>  7) | ((t[4] & 0x0f) << 4));
    r[11 * i + 6] = (unsigned char) ((t[4] >>  4) | ((t[5] & 0x01) << 7));
    r[11 * i + 7] = (unsigned char) ((t[5] >>  1) & 0xff);
    r[11 * i + 8] = (unsigned char) ((t[5] >>  9) | ((t[6] & 0x3f) << 2));
    r[11 * i + 9] = (unsigned char) ((t[6] >>  6) | ((t[7] & 0x07) << 5));
    r[11 * i + 10] = (unsigned char) ((t[7] >>  3));
  }

  for(j=0;j<NTRU_PACK_DEG-8*i;j++)
    t[j] = MODQ(a->coeffs[8*i+j]);
  for(; j<8; j++)
    t[j] = 0;

  switch(NTRU_PACK_DEG&0x07)
  {
    // cases 0 and 6 are impossible since 2 generates (Z/n)* and
    // p mod 8 in {1, 7} implies that 2 is a quadratic residue.
    case 4:
      r[11 * i + 0] = (unsigned char) (t[0]        & 0xff);
      r[11 * i + 1] = (unsigned char) (t[0] >>  8) | ((t[1] & 0x1f) << 3);
      r[11 * i + 2] = (unsigned char) (t[1] >>  5) | ((t[2] & 0x03) << 6);
      r[11 * i + 3] = (unsigned char) (t[2] >>  2) & 0xff;
      r[11 * i + 4] = (unsigned char) (t[2] >> 10) | ((t[3] & 0x7f) << 1);
      r[11 * i + 5] = (unsigned char) (t[3] >>  7) | ((t[4] & 0x0f) << 4);
      break;
    case 2:
      r[11 * i + 0] = (unsigned char) (t[0]        & 0xff);
      r[11 * i + 1] = (unsigned char) (t[0] >>  8) | ((t[1] & 0x1f) << 3);
      r[11 * i + 2] = (unsigned char) (t[1] >>  5) | ((t[2] & 0x03) << 6);
      break;
  }
}

void poly_Sq_frombytes(poly *r, const unsigned char *a)
{
  int i;
  for(i=0;i<NTRU_PACK_DEG/8;i++)
  {
    r->coeffs[8*i+0] = (a[11*i+ 0] >> 0) | (((uint16_t)a[11*i+ 1] & 0x07) << 8);
    r->coeffs[8*i+1] = (a[11*i+ 1] >> 3) | (((uint16_t)a[11*i+ 2] & 0x3f) << 5);
    r->coeffs[8*i+2] = (a[11*i+ 2] >> 6) | (((uint16_t)a[11*i+ 3] & 0xff) << 2) | (((uint16_t)a[11*i+ 4] & 0x01) << 10);
    r->coeffs[8*i+3] = (a[11*i+ 4] >> 1) | (((uint16_t)a[11*i+ 5] & 0x0f) << 7);
    r->coeffs[8*i+4] = (a[11*i+ 5] >> 4) | (((uint16_t)a[11*i+ 6] & 0x7f) << 4);
    r->coeffs[8*i+5] = (a[11*i+ 6] >> 7) | (((uint16_t)a[11*i+ 7] & 0xff) << 1) | (((uint16_t)a[11*i+ 8] & 0x03) <<  9);
    r->coeffs[8*i+6] = (a[11*i+ 8] >> 2) | (((uint16_t)a[11*i+ 9] & 0x1f) << 6);
    r->coeffs[8*i+7] = (a[11*i+ 9] >> 5) | (((uint16_t)a[11*i+10] & 0xff) << 3);
  }
  switch(NTRU_PACK_DEG&0x07)
  {
    // cases 0 and 6 are impossible since 2 generates (Z/n)* and
    // p mod 8 in {1, 7} implies that 2 is a quadratic residue.
    case 4:
      r->coeffs[8*i+0] = (a[11*i+ 0] >> 0) | (((uint16_t)a[11*i+ 1] & 0x07) << 8);
      r->coeffs[8*i+1] = (a[11*i+ 1] >> 3) | (((uint16_t)a[11*i+ 2] & 0x3f) << 5);
      r->coeffs[8*i+2] = (a[11*i+ 2] >> 6) | (((uint16_t)a[11*i+ 3] & 0xff) << 2) | (((uint16_t)a[11*i+ 4] & 0x01) << 10);
      r->coeffs[8*i+3] = (a[11*i+ 4] >> 1) | (((uint16_t)a[11*i+ 5] & 0x0f) << 7);
      break;
    case 2:
      r->coeffs[8*i+0] = (a[11*i+ 0] >> 0) | (((uint16_t)a[11*i+ 1] & 0x07) << 8);
      r->coeffs[8*i+1] = (a[11*i+ 1] >> 3) | (((uint16_t)a[11*i+ 2] & 0x3f) << 5);
      break;
  }
  r->coeffs[NTRU_N-1] = 0;
}
#endif

#ifdef PARAMS_3

void poly_Sq_tobytes(unsigned char *r, const poly *a)
{
  int i,j;
  uint16_t t[8];

  for(i=0;i<NTRU_PACK_DEG/8;i++)
  {
    for(j=0;j<8;j++)
      t[j] = MODQ(a->coeffs[8*i+j]);

    r[13 * i + 0] = (unsigned char) ( t[0]        & 0xff);
    r[13 * i + 1] = (unsigned char) ((t[0] >>  8) | ((t[1] & 0x07) << 5));
    r[13 * i + 2] = (unsigned char) ((t[1] >>  3) & 0xff);
    r[13 * i + 3] = (unsigned char) ((t[1] >> 11) | ((t[2] & 0x3f) << 2));
    r[13 * i + 4] = (unsigned char) ((t[2] >>  6) | ((t[3] & 0x01) << 7));
    r[13 * i + 5] = (unsigned char) ((t[3] >>  1) & 0xff);
    r[13 * i + 6] = (unsigned char) ((t[3] >>  9) | ((t[4] & 0x0f) << 4));
    r[13 * i + 7] = (unsigned char) ((t[4] >>  4) & 0xff);
    r[13 * i + 8] = (unsigned char) ((t[4] >> 12) | ((t[5] & 0x7f) << 1));
    r[13 * i + 9] = (unsigned char) ((t[5] >>  7) | ((t[6] & 0x03) << 6));
    r[13 * i + 10] = (unsigned char) ((t[6] >>  2) & 0xff);
    r[13 * i + 11] = (unsigned char) ((t[6] >> 10) | ((t[7] & 0x1f) << 3));
    r[13 * i + 12] = (unsigned char) ((t[7] >>  5));
  }

  for(j=0;j<NTRU_PACK_DEG-8*i;j++)
    t[j] = MODQ(a->coeffs[8*i+j]);
  for(; j<8; j++)
    t[j] = 0;

  switch(NTRU_PACK_DEG - 8*(NTRU_PACK_DEG/8))
  {
    // cases 0 and 6 are impossible since 2 generates (Z/n)* and
    // p mod 8 in {1, 7} implies that 2 is a quadratic residue.
    case 4:
      r[13 * i + 0] = (unsigned char) ( t[0]        & 0xff);
      r[13 * i + 1] = (unsigned char) (t[0] >>  8) | ((t[1] & 0x07) << 5);
      r[13 * i + 2] = (unsigned char) (t[1] >>  3) & 0xff;
      r[13 * i + 3] = (unsigned char) (t[1] >> 11) | ((t[2] & 0x3f) << 2);
      r[13 * i + 4] = (unsigned char) (t[2] >>  6) | ((t[3] & 0x01) << 7);
      r[13 * i + 5] = (unsigned char) (t[3] >>  1) & 0xff;
      r[13 * i + 6] = (unsigned char) (t[3] >>  9) | ((t[4] & 0x0f) << 4);
      break;
    case 2:
      r[13 * i + 0] = (unsigned char) ( t[0]        & 0xff);
      r[13 * i + 1] = (unsigned char) (t[0] >>  8) | ((t[1] & 0x07) << 5);
      r[13 * i + 2] = (unsigned char) (t[1] >>  3) & 0xff;
      r[13 * i + 3] = (unsigned char) (t[1] >> 11) | ((t[2] & 0x3f) << 2);
      break;
  }
}

void poly_Sq_frombytes(poly *r, const unsigned char *a)
{
  int i;
  for(i=0;i<NTRU_PACK_DEG/8;i++)
  {
    r->coeffs[8*i+0] =  a[13*i+ 0]       | (((uint16_t)a[13*i+ 1] & 0x1f) << 8);
    r->coeffs[8*i+1] = (a[13*i+ 1] >> 5) | (((uint16_t)a[13*i+ 2]       ) << 3) | (((uint16_t)a[13*i+ 3] & 0x03) << 11);
    r->coeffs[8*i+2] = (a[13*i+ 3] >> 2) | (((uint16_t)a[13*i+ 4] & 0x7f) << 6);
    r->coeffs[8*i+3] = (a[13*i+ 4] >> 7) | (((uint16_t)a[13*i+ 5]       ) << 1) | (((uint16_t)a[13*i+ 6] & 0x0f) <<  9);
    r->coeffs[8*i+4] = (a[13*i+ 6] >> 4) | (((uint16_t)a[13*i+ 7]       ) << 4) | (((uint16_t)a[13*i+ 8] & 0x01) << 12);
    r->coeffs[8*i+5] = (a[13*i+ 8] >> 1) | (((uint16_t)a[13*i+ 9] & 0x3f) << 7);
    r->coeffs[8*i+6] = (a[13*i+ 9] >> 6) | (((uint16_t)a[13*i+10]       ) << 2) | (((uint16_t)a[13*i+11] & 0x07) << 10);
    r->coeffs[8*i+7] = (a[13*i+11] >> 3) | (((uint16_t)a[13*i+12]       ) << 5);
  }
  switch(NTRU_PACK_DEG&0x07)
  {
    // cases 0 and 6 are impossible since 2 generates (Z/n)* and
    // p mod 8 in {1, 7} implies that 2 is a quadratic residue.
    case 4:
      r->coeffs[8*i+0] =  a[13*i+ 0]       | (((uint16_t)a[13*i+ 1] & 0x1f) << 8);
      r->coeffs[8*i+1] = (a[13*i+ 1] >> 5) | (((uint16_t)a[13*i+ 2]       ) << 3) | (((uint16_t)a[13*i+ 3] & 0x03) << 11);
      r->coeffs[8*i+2] = (a[13*i+ 3] >> 2) | (((uint16_t)a[13*i+ 4] & 0x7f) << 6);
      r->coeffs[8*i+3] = (a[13*i+ 4] >> 7) | (((uint16_t)a[13*i+ 5]       ) << 1) | (((uint16_t)a[13*i+ 6] & 0x0f) <<  9);
      break;
    case 2:
      r->coeffs[8*i+0] =  a[13*i+ 0]       | (((uint16_t)a[13*i+ 1] & 0x1f) << 8);
      r->coeffs[8*i+1] = (a[13*i+ 1] >> 5) | (((uint16_t)a[13*i+ 2]       ) << 3) | (((uint16_t)a[13*i+ 3] & 0x03) << 11);
      break;
  }
  r->coeffs[NTRU_N-1] = 0;
}

#endif

#ifdef PARAMS_4

void poly_Sq_tobytes(unsigned char *r, const poly *a)
{
  int i;

  for(i=0;i<NTRU_PACK_DEG/2;i++)
  {
    r[3 * i + 0] = (unsigned char) ( MODQ(a->coeffs[2 * i + 0]) & 0xff);
    r[3 * i + 1] = (unsigned char) ((MODQ(a->coeffs[2 * i + 0]) >>  8) | ((MODQ(a->coeffs[2 * i + 1]) & 0x0f) << 4));
    r[3 * i + 2] = (unsigned char) ((MODQ(a->coeffs[2 * i + 1]) >>  4));
  }
}

void poly_Sq_frombytes(poly *r, const unsigned char *a)
{
  int i;
  for(i=0;i<NTRU_PACK_DEG/2;i++)
  {
    r->coeffs[2*i+0] = (a[3*i+ 0] >> 0) | (((uint16_t)a[3*i+ 1] & 0x0f) << 8);
    r->coeffs[2*i+1] = (a[3*i+ 1] >> 4) | (((uint16_t)a[3*i+ 2] & 0xff) << 4);
  }
  r->coeffs[NTRU_N-1] = 0;
}

#endif
void poly_Rq_sum_zero_tobytes(unsigned char *r, const poly *a)
{
  poly_Sq_tobytes(r, a);
}

void poly_Rq_sum_zero_frombytes(poly *r, const unsigned char *a)
{
  int i;
  poly_Sq_frombytes(r,a);

  /* Set r[n-1] so that the sum of coefficients is zero mod q */
  r->coeffs[NTRU_N-1] = 0;
  for(i=0;i<NTRU_PACK_DEG;i++)
    r->coeffs[NTRU_N-1] -= r->coeffs[i];
}
