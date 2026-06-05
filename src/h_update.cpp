#include <stdio.h>
#include <math.h>
void h_update(double theta, double *wgt,
             int m1, double i_km1,double *zkm1,double *hkm1,
             int m2, double i_k,  double *zk,  double *hk
            )
{   double delta_k,x,rt_ik,rt_ikm1,rt_delta_k;
    int i,ii;
    delta_k=i_k-i_km1; /* incremental information */
    rt_delta_k=sqrt(delta_k);
    rt_ik=sqrt(i_k);
    rt_ikm1=sqrt(i_km1);
    for(i=0;i<=m2;i++)
    {   hk[i]=0.;
        for(ii=0;ii<=m1;ii++)
        {   x=(zk[i]*rt_ik-zkm1[ii]*rt_ikm1-theta*delta_k)/rt_delta_k;
            hk[i]+=hkm1[ii]*exp(-x*x/2)/2.506628275*rt_ik/rt_delta_k;
        }
        hk[i]*=wgt[i];
}   }
