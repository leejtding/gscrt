#define DEBUG 0
/* note: EXTREME_Z > 3 + log(r) +  Z(1-alpha) + Z(1-beta)
   per bottom of p 349 in Jennison and Turnbull */
#define EXTREME_Z 20
#define MAX_R 83
#include "R.h"
#include "Rmath.h"
#include "gs_design_crt.h"
/* Group sequential probability computation per Jennison & Turnbull
   x_n_anal- # of possible analyses in the group-sequential designs
            (interims + final)
   I     - statistical information available at each analysis
   a     - lower cutoff points for z statistic at each analysis (output)
   b     - upper cutoff points for z statistic at each analysis (output)
   prob_hi- input vector with probability of rejecting (Z>bj) at
           jth interim analysis, j=1...nanal
   prob_lo- input vector with probability of rejecting (Z<aj) at
           jth interim analysis, j=1...nanal
   x_tol  - relative change between iterations required to stop for 'convergence'
   x_r    - determinant of # of grid points for numerical integration
           r=17 will give a max of 201 points which is what they recommend
   retval- error flag returned; 0 if convergence; 1 indicates error
   printerr- 1 if error messages to be printed - other values suppress printing
*/
void gs_bound(int *x_n_anal,double *I,double *a,double *b,double *prob_lo,double *prob_hi,
             double *x_tol,int *x_r,int *retval,int *printerr) {
    int i,ii,j,m1,m2,r,n_anal;
    double plo,phi,dplo,dphi,btem=0.,atem=0.,atem2,btem2,rt_delta_k,rt_ik,rt_ikm1,xlo,xhi;
	double adelta,bdelta,tol;
    /* note: should allocat zwk & wwk dynamically...*/
    double zwk[1000],wwk[1000],hwk[1000],zwk2[1000],wwk2[1000],hwk2[1000],
           *z1,*z2,*w1,*w2,*h,*h2,*tem,rt_2pi;
    void h_1(double,int,double *,double,double *, double *);
    void h_update(double,double *,int,double,double *, double *,
                                 int,double,double *, double *);
    int grid_pts(int, double,double,double,double *, double *);
    r=x_r[0]; n_anal= x_n_anal[0]; tol=x_tol[0]; rt_2pi=2.506628274631;
    /* compute bounds at 1st interim analysis using inverse normal */
    if (n_anal < 1 || r<1 || r>MAX_R)
    {
        retval[0]=1;
 	 		if (*printerr)
			{	Rprintf("gs_bound error: illegal argument");
				if (n_anal<1) Rprintf("; n_anal=%d--must be > 0",n_anal);
				if (r<1 || r> MAX_R) Rprintf("; r=%d--must be >0 and <84",r);
				Rprintf("\n");
			}
	 		return;
	}
    if (prob_lo[0] <= 0) a[0] = -EXTREME_Z;
    else a[0]=qnorm(prob_lo[0],0.,1.,1,0);
    if (prob_hi[0] <= 0) b[0] = EXTREME_Z;
    else b[0]=qnorm(prob_hi[0],0.,1.,0,0);
/* set up work vectors */
    z1=zwk; w1=wwk; h=hwk;
    z2=zwk2; w2=wwk2; h2=hwk2;
    m1=grid_pts(r,0.,a[0],b[0],z1,w1);
    h_1(0.,m1,w1,I[0],z1,h);
    rt_ik=sqrt(I[0]);
    /* use Newton-Raphson to find subsequent interim analysis cutpoints */
    for(i=1;i<n_anal;i++)
    {   /* set up constants */
        rt_ikm1=rt_ik; rt_ik=sqrt(I[i]); rt_delta_k=sqrt(I[i]-I[i-1]);
        if (prob_lo[i]<=0.) atem2= -EXTREME_Z;
        else atem2=qnorm(prob_lo[i],0.,1.,1,0); 
        if (prob_hi[i]<=0.) btem2= EXTREME_Z;
        else btem2=qnorm(prob_hi[i],0.,1.,0,0);
        adelta=1.; bdelta=1.; j=0;
        while((adelta>tol || bdelta>tol) && j++ < EXTREME_Z)
	  {   plo=0.; phi=0.; dplo=0.; dphi=0.;
            atem=atem2; btem=btem2;
	/* compute probability of crossing boundaries & their derivatives */
            for(ii=0;ii<=m1;ii++)
            {   xlo=(z1[ii]*rt_ikm1-atem*rt_ik)/rt_delta_k;
                xhi=(z1[ii]*rt_ikm1-btem*rt_ik)/rt_delta_k;
                plo += h[ii]*pnorm(xlo,0.,1.,0,0);
                phi += h[ii]*pnorm(xhi,0.,1.,1,0);
                dplo+=h[ii]*exp(-xlo*xlo/2)/rt_2pi*rt_ik/rt_delta_k;
                dphi-=h[ii]*exp(-xhi*xhi/2)/rt_2pi*rt_ik/rt_delta_k;
            }
            /* use 1st order Taylor's series to update boundaries */
            /* maximum allowed change is 1 */
            /* maximum value allowed is z1[m1]*rt_ik to keep within grid points */
            adelta=prob_lo[i]-plo;
            if (adelta>dplo) atem2=atem+1.;
            else if (adelta < -dplo) atem2=atem-1.;
            else atem2=atem+(prob_lo[i]-plo)/dplo;
            if (atem2>EXTREME_Z) atem2=EXTREME_Z;
            else if (atem2 < -EXTREME_Z) atem2= -EXTREME_Z;
            bdelta=prob_hi[i]-phi;
            if (bdelta<dphi) btem2=btem+1.;
            else if (bdelta > -dphi) btem2=btem-1.;
            else btem2=btem+(prob_hi[i]-phi)/dphi;
            if (btem2>EXTREME_Z) btem2=EXTREME_Z;
            else if (btem2< -EXTREME_Z) btem2= -EXTREME_Z;
            if (atem2>btem2) atem2=btem2;
            adelta=atem2-atem; if (adelta<0) adelta= -adelta;
            bdelta=btem2-btem; if (bdelta<0) bdelta= -bdelta;
        }
        a[i]=atem; b[i]=btem;
/* if convergence did not occur, set flag for return value */
        if (adelta>tol ||bdelta > tol)
        {   if (*printerr) 
            {  Rprintf("gs_bound error: No convergence for boundary for interim %d; I=%7.0lf",i+1,I[i]);
	   			if (bdelta>tol) Rprintf("\n last 2 upper boundary values: %lf %lf\n",btem,btem2);
					if (adelta>tol) Rprintf("\n last 2 lower boundary values: %lf %lf\n",atem,atem);
				}
				retval[0]=1;
				return;
		  }
        if (i<n_anal-1)
        {   m2=grid_pts(r,0.,a[i],b[i],z2,w2);
            h_update(0.,w2,m1,I[i-1],z1,h,m2,I[i],z2,h2);
            m1=m2;
            tem=z1; z1=z2; z2=tem;
            tem=w1; w1=w2; w2=tem;
            tem=h;  h=h2;  h2=tem;
    }   }
    retval[0]=0; 
	 return;
}

