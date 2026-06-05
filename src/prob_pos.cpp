#include "R.h"
#include "Rmath.h"
double prob_pos(double theta,int m,double bk,double *z,double *h,double 
i_km1,double i_k)
{   int i;
    double xhi,prob,rt_ik,rt_ikm1,mu,rt_delta_k;
    mu=theta*(i_k-i_km1);
    rt_delta_k=sqrt(i_k-i_km1);
    rt_ik=sqrt(i_k);
    rt_ikm1=sqrt(i_km1);
    prob=0.;
    for(i=0;i<=m;i++)
    {   xhi=(z[i]*rt_ikm1+mu-bk*rt_ik)/rt_delta_k;
	    prob += pnorm(xhi,0.,1.,1,0)*h[i];
    }
    return(prob);
}

double prob_pos_2(double theta,int m,double bk,double *z,double *h,double 
i_km1,double i_k)
{   int i;
    double xhi1,xhi2,prob,rt_ik,rt_ikm1,mu,rt_delta_k;
    mu=theta*(i_k-i_km1);
    rt_delta_k=sqrt(i_k-i_km1);
    rt_ik=sqrt(i_k);
    rt_ikm1=sqrt(i_km1);
    prob=0.;
    for(i=0;i<=m;i++)
    {   xhi1=(z[i]*rt_ikm1+mu-bk*rt_ik)/rt_delta_k;
        xhi2=(z[i]*rt_ikm1+mu+bk*rt_ik)/rt_delta_k;
	    prob += (pnorm(xhi1,0.,1.,1,0)+pnorm(xhi2,0.,1.,0,0))*h[i];
    }
    return(prob);
}
