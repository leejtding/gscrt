#include "R.h"
#include "Rmath.h"
double prob_neg(double theta,int m,double ak,double *z,double *h,double 
i_km1,double i_k)
{   int i;
    double xlo,prob,rt_ik,rt_ikm1,mu,rt_delta_k;
    mu=theta*(i_k-i_km1);
    rt_delta_k=sqrt(i_k-i_km1);
    rt_ik=sqrt(i_k);
    rt_ikm1=sqrt(i_km1);
    prob=0.;
    for(i=0;i<=m;i++)
    {   xlo=(z[i]*rt_ikm1+mu-ak*rt_ik)/rt_delta_k;
	    prob += pnorm(xlo,0.,1.,0,0)*h[i];
    }
    return(prob);
}

double prob_neg_2(double theta,int m,double ak,double *z,double *h,double 
i_km1,double i_k)
{   int i;
    double xlo1,xlo2,prob,rt_ik,rt_ikm1,mu,rt_delta_k;
    mu=theta*(i_k-i_km1);
    rt_delta_k=sqrt(i_k-i_km1);
    rt_ik=sqrt(i_k);
    rt_ikm1=sqrt(i_km1);
    prob=0.;
    for(i=0;i<=m;i++)
    {   xlo1=(z[i]*rt_ikm1+mu-ak*rt_ik)/rt_delta_k;
        xlo2=(z[i]*rt_ikm1+mu+ak*rt_ik)/rt_delta_k;
	    prob += (pnorm(xlo1,0.,1.,0,0)-pnorm(xlo2,0.,1.,0,0))*h[i];
    }
    return(prob);
}
