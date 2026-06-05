#define DEBUG 0
#define EXTREME_Z 20
#define MAX_R 83
#include "R.h"
#include "Rmath.h"
#include "gs_design_crt.h"
/* Group sequential probability computation per Jennison & Turnbull
   Computes upper bound to have input crossing probabilities given fixed input lower bound.
   x_n_anal- # of possible analyses in the group-sequential designs
           (interims + final)
   x_theta- drift parameter
   I     - statistical information available at each analysis
   a     - lower cutoff points for z statistic at each analysis (input)
   b     - upper cutoff points for z statistic at each analysis (output)
   prob_lo- output vector with probability of rejecting (Z<aj) at
           jth interim analysis, j=1...nanal
   prob_hi- input vector with probability of rejecting (Z>bj) at
           jth interim analysis, j=1...nanal
   tol   - relative change between iterations required to stop for 'convergence'
   x_r    - controls # of grid points for numerical integration per Jennison & Turnbull
	        they recommend r=17 (r=18 is default - a little more accuracy)
   retval- error flag returned; 0 if convergence; 1 indicates error
   printerr- 1 if error messages to be printed - other values suppress printing
*/
void gs_upper_1(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
              double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr) {
    int i,ii,j,m1,m2,r,n_anal;
    double plo=0.,phi,dphi,btem=0.,btem2,rt_delta_k,rt_ik,rt_ikm1,xlo,xhi,theta,mu,tol,bdelta;
    /* note: should allocat zwk & wwk dynamically...*/
    double zwk[1000],wwk[1000],hwk[1000],zwk2[1000],wwk2[1000],hwk2[1000],
           *z1,*z2,*w1,*w2,*h,*h2,*tem,rt_2pi;
    void h_1(double,int,double *,double,double *, double *);
    void h_update(double,double *,int,double,double *, double *,
                                 int,double,double *, double *);
    int grid_pts(int,double,double,double,double *, double *);
    r=x_r[0]; n_anal= x_n_anal[0]; theta= x_theta[0]; tol=x_tol[0]; rt_2pi=2.506628274631;
    if (n_anal<1 || r<1 || r>MAX_R) 
	{	retval[0]=1;
        if (*printerr)
        {	Rprintf("gs_upper error: illegal argument");
            if (n_anal<1) Rprintf("; n_anal=%d--must be > 0",n_anal);
            if (r<1 || r> MAX_R) Rprintf("; r=%d--must be >0 and <84",r);
            Rprintf("\n");
        }
        return;
	}
    rt_ik=sqrt(I[0]); mu=rt_ik*theta;
    prob_lo[0]=pnorm(mu-a[0],0.,1.,0,0);
    if (prob_hi[0] <= 0.) b[0]=EXTREME_Z;
    else b[0]=qnorm(prob_hi[0],mu,1.,0,0);
    if (n_anal==1) {retval[0]=0; return;}
    
    /* set up work vectors */
    z1=zwk; w1=wwk; h=hwk;
    z2=zwk2; w2=wwk2; h2=hwk2;
	if (DEBUG) Rprintf("r=%d mu=%lf a[0]=%lf b[0]=%lf\n",r,mu,a[0],b[0]);
    m1=grid_pts(r,mu,a[0],b[0],z1,w1);
    h_1(theta,m1,w1,I[0],z1,h);

    /* use Newton-Raphson to find subsequent interim analysis cutpoints */
	retval[0]=0;
    for(i=1;i<n_anal;i++)
    {   rt_ikm1=rt_ik; rt_ik=sqrt(I[i]); mu=rt_ik*theta; rt_delta_k=sqrt(I[i]-I[i-1]);
        if (prob_hi[i] <= 0.) btem2=EXTREME_Z;
        else btem2=qnorm(prob_hi[i],mu,1.,0,0); 
        bdelta=1.; j=0;
        while((bdelta>tol) && j++ < 20)
		{   phi=0.; dphi=0.; plo=0.;
            btem=btem2;
			if (DEBUG) Rprintf("i=%d m1=%d\n",i,m1);
	        
            /* compute probability of crossing boundaries & their derivatives */
            for(ii=0;ii<=m1;ii++)
            {   xhi=(z1[ii]*rt_ikm1-btem*rt_ik+theta*(I[i]-I[i-1]))/rt_delta_k;
                phi += pnorm(xhi,0.,1.,1,0)*h[ii];
				xlo=(z1[ii]*rt_ikm1-a[i]*rt_ik+theta*(I[i]-I[i-1]))/rt_delta_k;
				plo += pnorm(xlo,0.,1.,0,0)*h[ii];
                dphi-=h[ii]*exp(-xhi*xhi/2)/rt_2pi*rt_ik/rt_delta_k;
				if (DEBUG) Rprintf("m1=%d ii=%d xhi=%lf phi=%lf xlo=%lf plo=%lf dphi=%lf\n",m1,ii,xhi,phi,xlo,plo,dphi);
            }

            /* use 1st order Taylor's series to update boundaries */
            /* maximum allowed change is 1 */
            /* maximum value allowed is EXTREME_Z */
            if (DEBUG) Rprintf("i=%2d j=%2d plo=%lf btem=%lf phi=%lf dphi=%lf\n",i,j,plo,btem,phi,dphi);
            bdelta=prob_hi[i]-phi;
            if (bdelta<dphi) btem2=btem+1.;
            else if (bdelta > -dphi) btem2=btem-1.;
            else btem2=btem+(prob_hi[i]-phi)/dphi;
            if (btem2>EXTREME_Z) btem2=EXTREME_Z;
            else if (btem2< -EXTREME_Z) btem2= -EXTREME_Z;

            /* adjust boundaries as needed */
            bdelta=btem2-btem; if (bdelta<0) bdelta= -bdelta;
        }
        b[i]=btem;
        prob_lo[i]=plo;
	    
        /* if convergence did not occur, set flag for return value */
        if (bdelta > tol)
		{   if (*printerr) Rprintf("gs_upper error: No convergence for boundary for interim %d; I=%7.0lf; last 2 upper boundary values: %lf %lf\n",
				i+1,I[i],btem,btem2);
			retval[0]=1;
            return;
		}
        if (i<n_anal-1)
        {   m2=grid_pts(r,mu,a[i],b[i],z2,w2);
            h_update(theta,w2,m1,I[i-1],z1,h,m2,I[i],z2,h2);
            m1=m2;
            tem=z1; z1=z2; z2=tem;
            tem=w1; w1=w2; w2=tem;
            tem=h;  h=h2;  h2=tem;
        }
    }
    return;
}
