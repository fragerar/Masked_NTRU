#include "masking.h"
#include "random.h"

static uint32_t genrand(int l)
{
  if (l==32) return rand32();
  return rand32() & ((1 << l)-1);
}



static void Expand(uint32_t *x,uint32_t *xp,int k,int n2,int n)
{
  for(int i=0;i<n/2;i++)
  {
    uint32_t r=genrand(k);
    xp[2*i]=x[i] ^ r;
    xp[2*i+1]=r;
  }
  if ((n & 1)==1) 
  {
    if (n2==n/2)
      xp[n-1]=0;
    else
      xp[n-1]=x[n2-1];
  }
}

static void SecAnd(uint32_t *a,uint32_t *b,uint32_t *c,int k,int n)
{
  for(int i=0;i<n;i++)
    c[i]=a[i] & b[i];

  for(int i=0;i<n;i++)
  {
    for(int j=i+1;j<n;j++)
    {
      uint32_t tmp=rand32(); //rand();
      uint32_t tmp2=(tmp ^ (a[i] & b[j])) ^ (a[j] & b[i]);
      c[i]^=tmp;
      c[j]^=tmp2;
    }
  }
  for(int i=0;i<n;i++) c[i]=c[i] % (1 << k);
}

static void SecAdd(uint32_t *x,uint32_t *y,uint32_t *z,int k,int n)
{
  uint32_t u[n];
  for(int i=0;i<n;i++) u[i]=0;
  uint32_t w[n];
  SecAnd(x,y,w,k,n);
  uint32_t a[n];
  for(int i=0;i<n;i++) a[i]=x[i] ^ y[i];
  for(int j=0;j<k-1;j++)
  {
    uint32_t ua[n];
    SecAnd(u,a,ua,k,n);
    for(int i=0;i<n;i++) u[i]=(2*(ua[i] ^ w[i])) % (1 << k);
  }
  for(int i=0;i<n;i++) z[i]=x[i] ^ y[i] ^ u[i];
}


static uint32_t GoubinAB(uint32_t A,uint32_t r,int k)
{
  uint32_t G=rand32();
  uint32_t T=G << 1;
  uint32_t x=G ^ r;
  uint32_t O=G & x;
  x=T ^ A;
  G=G ^ x;
  G=G & r;
  O=O ^ G;
  G=T & A;
  O=O ^ G;
  for(int i=1;i<k;i++)
  {
    G=T & r;
    G=G ^ O;
    T=T & A;
    G=G ^ T;
    T=G << 1;
  }
  x=x ^ T;
  return x;
}


static void ConvertAB(uint32_t *A,uint32_t *z,int k,int n)
{
  if(n==1)
  {
    z[0]=A[0];
    return;
  }

  if(n==2)
  {
    z[0]=GoubinAB(A[0],A[1],k);
    z[1]=A[1];
    return;
  }

  uint32_t x[n/2];
  ConvertAB(A,x,k,n/2);
  uint32_t xp[n];
  Expand(x,xp,k,n/2,n);
  
  uint32_t y[(n+1)/2];
  ConvertAB(A+n/2,y,k,(n+1)/2);
  uint32_t yp[n];
  Expand(y,yp,k,(n+1)/2,n);

  SecAdd(xp,yp,z,k,n);
}


void convert_A2B_CGV14(Masked* y, const Masked* x, unsigned k){
  uint32_t a[N_SHARES], b[N_SHARES];
  for(int i=0; i < N_SHARES; ++i) a[i] = (uint32_t) x->shares[i];
  ConvertAB(a, b, k, N_SHARES);
  for(int i=0; i < N_SHARES; ++i) y->shares[i] = (uint16_t) b[i]%(1<<k);
}

static void print_CGV(uint32_t* y){
  int t=0;
  printf(" (");
  for(int i=0; i < MASKING_ORDER; ++i){
    printf("%u, ",y[i]);
    t ^= y[i];
  }
  t ^= y[MASKING_ORDER];
  printf("%i) = (", y[MASKING_ORDER]);
  for(int i=0; i < MASKING_ORDER; ++i){
    printf("0x%X, ", y[i]);
  }
  printf("0x%x) = %u\n", y[MASKING_ORDER], t);
}

