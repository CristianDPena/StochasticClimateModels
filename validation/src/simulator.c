#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
/* Coupled Euler paths: every observation uses the same sixteen fine Brownian
   increments; levels aggregate these to 4, 8, or 16 latent substeps. */
SEXP coupled_latent(SEXP nn, SEXP dd, SEXP rr, SEXP ll) {
  int n=asInteger(nn), L=LENGTH(ll); double dt=asReal(dd), rho=asReal(rr);
  if(n<1 || dt<=0 || rho<=0) error("Invalid simulator inputs");
  SEXP out=PROTECT(allocMatrix(REALSXP,n+1,L));
  double *v=REAL(out), state[3]={0,0,0};
  if(L<1 || L>3) error("Invalid levels");
  for(int l=0;l<L;l++){int m=INTEGER(ll)[l];if(m!=4&&m!=8&&m!=16)error("Invalid level");v[l*(R_xlen_t)(n+1)]=0;}
  GetRNGstate();
  for(int i=0;i<n;i++) {
    double z[16]; for(int j=0;j<16;j++) z[j]=norm_rand()*sqrt(dt/16.0);
    for(int l=0;l<L;l++) {
      int m=INTEGER(ll)[l], d=16/m;
      for(int j=0;j<m;j++){double dw=0;for(int a=0;a<d;a++)dw+=z[j*d+a];
        state[l] += -rho*rho*state[l]*dt/m + rho*sqrt(1+state[l]*state[l])*dw;}
      v[(R_xlen_t)l*(n+1)+i+1]=state[l];
    }
    if(i%100000==0) R_CheckUserInterrupt();
  }
  PutRNGstate(); UNPROTECT(1); return out;
}
