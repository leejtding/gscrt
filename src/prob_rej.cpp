#include "R.h"
#include "Rmath.h"
#include "gs_design_crt.h"
/* Group sequential boundary crossing probability computation per Jennison & Turnbull
   This version uses all pointer arguments so that it can be called from R or Splus
   x_n_anal - # of possible analyses in the group-sequential designs
            (interims + final)
	n_theta - # of theta values for which boundary crossing probabilities are to be computed
   theta  - vector of drift parameters
   I      - statistical information available at each analysis
   a      - lower cutoff points for z statistic at each analysis
   b      - upper cutoff points for z statistic at each analysis
   x_prob_hi- vector to return probability of rejecting (Z>bj) at
            jth interim analysis, j=1...nanal
   x_prob_lo- vector to return probability of rejecting (Z<aj) at
            jth interim analysis, j=1...nanal
   x_r     - determinant of # of grid points for numerical integration
            r=17 will give a max of 201 points which is what they recommend
*/
void prob_rej_1(int *x_n_anal,int *n_theta,double *x_theta,double *I,double *a,double *b,
              double *x_prob_lo,double *x_prob_hi,int *x_r)
{   int r,i,m1,m2,n_anal,k;
    double theta;
    double *prob_lo,*prob_hi;
    double prob_neg(double,int,double,double *,double *,double,double);
    double prob_pos(double,int,double,double *,double *,double,double);
/* note: should allocat zwk & wwk dynamically...*/
    double mu,zwk[1000],wwk[1000],hwk[1000],zwk2[1000],wwk2[1000],hwk2[1000],
           *z1,*z2,*w1,*w2,*h,*h2,*tem;
    void h_1(double,int,double *,double,double *, double *);
    void h_update(double,double *,int,double,double *, double *,
                                 int,double,double *, double *);
    int grid_pts(int, double,double,double,double *, double *);
    r=x_r[0]; n_anal=x_n_anal[0]; 
    for(k=0;k<n_theta[0];k++)
    {  theta=x_theta[k];
       prob_lo=x_prob_lo+k*n_anal;
       prob_hi=x_prob_hi+k*n_anal;
/* compute probability of rejecting at 1st interim analysis */
       if (n_anal < 1) return;
       mu=theta*sqrt(I[0]);
       prob_lo[0]=pnorm(a[0],mu,1.,1,0);
       prob_hi[0]=pnorm(b[0],mu,1.,0,0);
/* compute h_1 */
       z1=zwk; w1=wwk; h=hwk;
       m1=grid_pts(r,mu,a[0],b[0],z1,w1);
       h_1(theta,m1,w1,I[0],z1,h);
       z2=zwk2; w2=wwk2; h2=hwk2;
/* update h and compute rejection probabilities at each interim */
       for(i=1;i<n_anal;i++)
       {   prob_hi[i]=prob_pos(theta,m1,b[i],z1,h,I[i-1],I[i]);
           prob_lo[i]=prob_neg(theta,m1,a[i],z1,h,I[i-1],I[i]);
           if (i<n_anal-1)
	   {   mu=theta*sqrt(I[i]);
	       m2=grid_pts(r,mu,a[i],b[i],z2,w2);
              h_update(theta,w2,m1,I[i-1],z1,h,m2,I[i],z2,h2);
              m1=m2;
              tem=z1; z1=z2; z2=tem;
              tem=w1; w1=w2; w2=tem;
              tem=h;  h=h2;  h2=tem;
}   }   }  }

void prob_rej_2(int *x_n_anal,int *n_theta,double *x_theta,double *I,double *a,double *b,
              double *x_prob_lo,double *x_prob_hi,int *x_r)
{   int r,i,m1,m2,n_anal,k;
    double theta;
    double *prob_lo,*prob_hi;
    double prob_neg_2(double,int,double,double *,double *,double,double);
    double prob_pos_2(double,int,double,double *,double *,double,double);
/* note: should allocat zwk & wwk dynamically...*/
    double mu,zwk[1000],wwk[1000],hwk[1000],zwk2[1000],wwk2[1000],hwk2[1000],
           *z1,*z2,*w1,*w2,*h,*h2,*tem;
    void h_1(double,int,double *,double,double *, double *);
    void h_update(double,double *,int,double,double *, double *,
                                 int,double,double *, double *);
    int grid_pts_2(int, double,double,double,double *, double *);
    r=x_r[0]; n_anal=x_n_anal[0]; 
    for(k=0;k<n_theta[0];k++)
    {  theta=x_theta[k];
       prob_lo=x_prob_lo+k*n_anal;
       prob_hi=x_prob_hi+k*n_anal;
/* compute probability of rejecting at 1st interim analysis */
       if (n_anal < 1) return;
       mu=theta*sqrt(I[0]);
       prob_lo[0]=pnorm(a[0],mu,1.,1,0)-pnorm(-a[0],mu,1.,1,0);
       prob_hi[0]=pnorm(b[0],mu,1.,0,0)+pnorm(-b[0],mu,1.,1,0);
/* compute h_1 */
       z1=zwk; w1=wwk; h=hwk;
       m1=grid_pts_2(r,mu,a[0],b[0],z1,w1);
       h_1(theta,m1,w1,I[0],z1,h);
       z2=zwk2; w2=wwk2; h2=hwk2;
/* update h and compute rejection probabilities at each interim */
       for(i=1;i<n_anal;i++)
       {   prob_hi[i]=prob_pos_2(theta,m1,b[i],z1,h,I[i-1],I[i]);
           prob_lo[i]=prob_neg_2(theta,m1,a[i],z1,h,I[i-1],I[i]);
           if (i<n_anal-1)
	   {   mu=theta*sqrt(I[i]);
	       m2=grid_pts_2(r,mu,a[i],b[i],z2,w2);
              h_update(theta,w2,m1,I[i-1],z1,h,m2,I[i],z2,h2);
              m1=m2;
              tem=z1; z1=z2; z2=tem;
              tem=w1; w1=w2; w2=tem;
              tem=h;  h=h2;  h2=tem;
}   }   }  }
