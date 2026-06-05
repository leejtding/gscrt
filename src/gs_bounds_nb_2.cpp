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
   prob_lo- input vector with probability of rejecting (Z<aj) at
           jth interim analysis, j=1...nanal
   prob_hi- input vector with probability of rejecting (Z>bj) at
           jth interim analysis, j=1...nanal
   x_tol  - relative change between iterations required to stop for 'convergence'
   x_r    - determinant of # of grid points for numerical integration
           r=17 will give a max of 201 points which is what they recommend
   retval- error flag returned; 0 if convergence; 1 indicates error
   printerr- 1 if error messages to be printed - other values suppress printing
*/
void gs_bounds_nb_2(int *x_n_anal, double *x_theta, double *I, double *a, double *b,
                 double *prob_lo, double *prob_hi, double *x_tol, int *x_r, int *retval,
                 int *printerr) {
    int i,i1,i2,j,m11,m12,m21,m22,r,n_anal;
    double plo,phi,dplo,dphi,btem=0.,atem=0.,atem2,btem2,rt_delta_k,rt_ik,rt_ikm1,xlo,xlo2,xhi,theta,mu1,mu2;
	double adelta,bdelta,tol;
    /* note: should allocat zwk & wwk dynamically...*/
    double zwk11[1000],wwk11[1000],hwk11[1000],zwk12[1000],wwk12[1000],hwk12[1000],
           zwk21[1000],wwk21[1000],hwk21[1000],zwk22[1000],wwk22[1000],hwk22[1000],
           *z11,*z12,*w11,*w12,*h11,*h12,*z21,*z22,*w21,*w22,*h21,*h22,*tem,rt_2pi;
    void h_1(double,int,double *,double,double *,double *);
    void h_update(double,double *,int,double,double *,double *,int,double,double *,double *);
    int grid_pts_2(int,double,double,double,double *,double *);
    r=x_r[0]; n_anal=x_n_anal[0]; theta=x_theta[0]; tol=x_tol[0]; rt_2pi=2.506628274631;

    /* compute bounds at 1st interim analysis using inverse normal */
    if (n_anal<1 || r<1 || r>MAX_R)
    {   retval[0]=1;
        if (*printerr)
        {	Rprintf("gs_bound error: illegal argument");
            if (n_anal<1) Rprintf("; n_anal=%d--must be > 0",n_anal);
            if (r<1 || r> MAX_R) Rprintf("; r=%d--must be >0 and <84",r);
            Rprintf("\n");
        }
        return;
	}
    rt_ik=sqrt(I[0]); mu1=0.; mu2=rt_ik*theta;
    if (prob_hi[0] <= 0.) b[0] = EXTREME_Z;
    else b[0]=qnorm(prob_hi[0],mu1,1.,0,0);
    if (prob_lo[0]<=0.) atem2= 0.;
    else atem2=qnorm(1-prob_lo[0],mu2,1.,0,0);
    adelta=1.; j=0;
    while((adelta>tol) && j++ < EXTREME_Z)
    {   plo=0.; dplo=0.; atem=atem2;

        /* construct lower boundary */
        /* compute probability of crossing lower boundaries & their derivatives under H_1 */
        xlo=atem-mu2;
        xlo2=-atem-mu2;
        plo+=(pnorm(xlo,0.,1.,1,0)-pnorm(xlo2,0.,1.,1,0));
        dplo+=(exp(-xlo*xlo/2)+exp(-xlo2*xlo2/2))/rt_2pi;

        /* use 1st order Taylor's series to update upper boundary under H_1 */
        /* maximum allowed change is 1 */
        /* maximum value allowed is z1[m1]*rt_ik to keep within grid points */
        adelta=prob_lo[0]-plo;
        if (adelta>dplo) atem2=atem+1.;
        else if (adelta < -dplo) atem2=atem-1.;
        else atem2=atem+(prob_lo[0]-plo)/dplo;
        if (atem2>EXTREME_Z) atem2=EXTREME_Z;
        else if (atem2 < 0) atem2= 0;
        adelta=atem2-atem; if (adelta<0) adelta= -adelta;
    }
    a[0]=atem;
    if (n_anal==1) {retval[0]=0; return;}

    /* set up work vectors */
    z11=zwk11; w11=wwk11; h11=hwk11; z12=zwk12; w12=wwk12; h12=hwk12;
    z21=zwk21; w21=wwk21; h21=hwk21; z22=zwk22; w22=wwk22; h22=hwk22;
    
    m11=grid_pts_2(r,mu1,0,b[0],z11,w11);
    h_1(0.,m11,w11,I[0],z11,h11);
    
    m12=grid_pts_2(r,mu2,a[0],b[0],z12,w12);
    h_1(theta,m12,w12,I[0],z12,h12);

    /* use Newton-Raphson to find subsequent interim analysis cutpoints */
    if (*printerr) Rprintf("Start: r=%d mu1=%lf mu2=%lf a[0]=%lf b[0]=%lf\n",r,mu1,mu2,a[0],b[0]);
    retval[0]=0;
    for(i=1;i<n_anal;i++)
    {   /* set up constants */
        rt_ikm1=rt_ik; rt_ik=sqrt(I[i]); mu2=rt_ik*theta; rt_delta_k=sqrt(I[i]-I[i-1]);
        if (prob_lo[i]<=0.) atem2= 0.;
        else atem2=qnorm(1-prob_lo[i],mu2,1.,0,0);
        if (prob_hi[i]<=0.) btem2= EXTREME_Z;
        else btem2=qnorm(prob_hi[i],mu1,1.,0,0);
        
        /* find upper boundary */
        bdelta=1.; j=0;
        if (*printerr) Rprintf("desired prob_hi=%lf\n", prob_hi[i]);
        while((bdelta>tol) && j++ < EXTREME_Z)
	    {   phi=0.; dphi=0.; btem=btem2;
            if (*printerr) Rprintf("i=%d m11=%d m12=%d\n",i,m11,m12);
	        
            /* construct upper boundary */
            /* compute probability of crossing upper boundaries & their derivatives under H_0 */
            for(i1=0;i1<=m11;i1++)
            {   xhi=(z11[i1]*rt_ikm1-btem*rt_ik)/rt_delta_k;
                phi+=h11[i1]*pnorm(xhi,0.,1.,1,0);
                dphi-=h11[i1]*exp(-xhi*xhi/2)/rt_2pi*rt_ik/rt_delta_k;
            }

            /* use 1st order Taylor's series to update upper boundary under H_0*/
            /* maximum allowed change is 1 */
            /* maximum value allowed is z1[m1]*rt_ik to keep within grid points */
            if (*printerr) Rprintf("i=%2d j=%2d btem=%lf phi=%lf dphi=%lf\n",i,j,btem,phi,dphi);       
            bdelta=prob_hi[i]-phi;
            if (bdelta<dphi) btem2=btem+1.;
            else if (bdelta > -dphi) btem2=btem-1.;
            else btem2=btem+(prob_hi[i]-phi)/dphi;
            if (btem2>EXTREME_Z) btem2=EXTREME_Z;
            else if (btem2< -EXTREME_Z) btem2= -EXTREME_Z;
            bdelta=btem2-btem; if (bdelta<0) bdelta= -bdelta;
        }
        b[i]=btem;
        
        /* find lower boundary */
        adelta=1.; j=0;
        if (*printerr) Rprintf("desired prob_lo=%lf\n", prob_lo[i]);
        while((adelta>tol) && j++ < EXTREME_Z)
	    {   plo=0.; dplo=0.; atem=atem2;
            if (*printerr) Rprintf("i=%d m11=%d m12=%d\n",i,m11,m12);

            /* construct lower boundary */
            /* compute probability of crossing lower boundaries & their derivatives under H_1 */
            for(i2=0;i2<=m12;i2++)
            {   
                xlo=(-z12[i2]*rt_ikm1+atem*rt_ik-theta*(I[i]-I[i-1]))/rt_delta_k;
                xlo2=(-z12[i2]*rt_ikm1-atem*rt_ik-theta*(I[i]-I[i-1]))/rt_delta_k;
                plo+=h12[i2]*(pnorm(xlo,0.,1.,1,0)-pnorm(xlo2,0.,1.,1,0));
                dplo+=h12[i2]*(rt_ik/rt_delta_k)*(exp(-xlo*xlo/2)+exp(-xlo2*xlo2/2))/rt_2pi;
            }

            /* use 1st order Taylor's series to update upper boundary under H_1 */
            /* maximum allowed change is 1 */
            /* maximum value allowed is z1[m1]*rt_ik to keep within grid points */
            if (*printerr) Rprintf("i=%2d j=%2d atem=%lf plo=%lf dplo=%lf\n",i,j,atem,plo,dplo);
            adelta=prob_lo[i]-plo;
            if (adelta>dplo) atem2=atem+1.;
            else if (adelta < -dplo) atem2=atem-1.;
            else atem2=atem+(prob_lo[i]-plo)/dplo;
            if (atem2>EXTREME_Z) atem2=EXTREME_Z;
            adelta=atem2-atem; if (adelta<0) adelta= -adelta;
        }
        a[i]=atem;
        
        /* if convergence did not occur, set flag for return value */
        if (adelta>tol || bdelta > tol)
        {   if (*printerr) 
            {  Rprintf("gs_bound error: No convergence for boundary for interim %d; I=%7.0lf",i+1,I[i]);
	   			if (bdelta>tol) Rprintf("\n last 2 upper boundary values: %lf %lf\n",btem,btem2);
				if (adelta>tol) Rprintf("\n last 2 lower boundary values: %lf %lf\n",atem,atem2);
			}
			retval[0]=1;
		}
        if (i<n_anal-1)
        {   m21=grid_pts_2(r,mu1,0,b[i],z21,w21);
            m22=grid_pts_2(r,mu2,a[i],b[i],z22,w22);
            h_update(0.,w21,m11,I[i-1],z11,h11,m21,I[i],z21,h21);
            h_update(theta,w22,m12,I[i-1],z12,h12,m22,I[i],z22,h22);
            m11=m21; m12=m22;
            tem=z11; z11=z21; z21=tem;
            tem=w11; w11=w21; w21=tem;
            tem=h11; h11=h21; h21=tem;
            tem=z12; z12=z22; z22=tem;
            tem=w12; w12=w22; w22=tem;
            tem=h12; h12=h22; h22=tem;
        }
    }
    retval[0]=0; 
	return;
}

