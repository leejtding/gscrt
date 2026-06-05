#include "R.h"
#include "Rmath.h"
#include "gs_design_crt.h"
/* sub-density function (integrates to < 1) between bounds at interim
   analysis for group sequential design */
void gs_density(double *den, int *x_n_anal, int *n_theta, double *x_theta,
			   double *I, double *a, double *b, double *x_z,
			   int *z_len, int *x_r)
{   int r, i, j, k, m1, m2, n_anal, nz;
    double z, mu, theta;
    double *zwk,*wwk,*hwk,*zwk2,*wwk2,*hwk2;
    double *z1,*z2,*w1,*w2,*h,*h2,*tem;
    void h_1(double,int,double *,double,double *, double *);
    void h_update(double,double *,int,double,double *, double *,
                                 int,double,double *, double *);
    int grid_pts(int, double,double,double,double *, double *);
    r=x_r[0]; n_anal=x_n_anal[0]; nz=z_len[0];
/* if density is at 1st analysis, just return normal density */
    if (n_anal < 1) return;
    if (n_anal == 1)
    {   j = 0;
        for(k=0; k < n_theta[0]; k++)
        {   mu = sqrt(I[0]) * x_theta[k];
            for(i=0; i < nz; i++) 
            {   z = x_z[i] - mu;
                den[j++] = exp(-z*z/2)/2.506628275;
        }   }
        return;
    } 
/* otherwise, compute density like in prob_rej */
	zwk = (double *) R_alloc(r * 12 - 3, sizeof(double));
	wwk = (double *) R_alloc(r * 12 - 3, sizeof(double));
	hwk = (double *) R_alloc(r * 12 - 3, sizeof(double));
	zwk2 = (double *) R_alloc(r * 12 - 3, sizeof(double));
	wwk2 = (double *) R_alloc(r * 12 - 3, sizeof(double));
	hwk2 = (double *) R_alloc(r * 12 - 3, sizeof(double));
    for(k=0; k < n_theta[0]; k++)
    {   theta=x_theta[k];
        mu = theta * sqrt(I[0]);
/* compute h_1 */
        z1=zwk; w1=wwk; h=hwk;
        m1=grid_pts(r,mu,a[0],b[0],z1,w1);
        h_1(theta,m1,w1,I[0],z1,h);
        z2=zwk2; w2=wwk2; h2=hwk2;
/* update h and compute rejection probabilities at each interim */
       for(i=1;i<n_anal;i++)
       {   mu=theta*sqrt(I[i]);
           if (i < n_anal - 1) 
			   m2 = grid_pts(r, mu, a[i], b[i], z2, w2);
           else
           {   m2 = nz - 1; z2 = x_z; h2 = den + k * nz;
               for(j = 0; j < nz; j++) w2[j]=1.;
           }
           h_update(theta,w2,m1,I[i-1],z1,h,m2,I[i],z2,h2);
           m1=m2;
           tem=z1; z1=z2; z2=tem;
           tem=w1; w1=w2; w2=tem;
           tem=h;  h=h2;  h2=tem;
       }
    }
}
