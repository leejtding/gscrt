#ifndef GS_DESIGN_CRT_H
#define GS_DESIGN_CRT_H

#ifdef __cplusplus
extern "C" {
#endif

void gs_bound(int *x_n_anal,double *I,double *a,double *b,double *prob_lo,double *prob_hi,
              double *x_tol,int *x_r,int *retval,int *printerr);
void gs_upper_1(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
              double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr);
void gs_upper_2(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
              double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr);
void gs_lower_1(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
              double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr);
void gs_lower_2(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
              double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr);
void gs_bounds_1(int *x_n_anal, double *x_theta, double *I, double *a, double *b,
               double *prob_lo, double *prob_hi, double *x_tol, int *x_r, int *retval,
               int *printerr);
void gs_bounds_2(int *x_n_anal, double *x_theta, double *I, double *a, double *b,
               double *prob_lo, double *prob_hi, double *x_tol, int *x_r, int *retval,
               int *printerr);
void gs_bounds_nb_1(int *x_n_anal, double *x_theta, double *I, double *a, double *b,
                 double *prob_lo, double *prob_hi, double *x_tol, int *x_r, int *retval,
                 int *printerr);
void gs_bounds_nb_2(int *x_n_anal, double *x_theta, double *I, double *a, double *b,
                 double *prob_lo, double *prob_hi, double *x_tol, int *x_r, int *retval,
                 int *printerr);
void prob_rej_1(int *x_n_anal,int *n_theta,double *x_theta,double *I,double *a,double *b,
              double *x_prob_lo,double *x_prob_hi,int *x_r);
void prob_rej_2(int *x_n_anal,int *n_theta,double *x_theta,double *I,double *a,double *b,
              double *x_prob_lo,double *x_prob_hi,int *x_r);
void prob_stop_1(int *x_n_anal,int *n_theta,double *x_theta,double *I,double *a,double *b,
               double *x_prob_lo,double *x_prob_hi,int *x_r);
void prob_stop_2(int *x_n_anal,int *n_theta,double *x_theta,double *I,double *a,double *b,
               double *x_prob_lo,double *x_prob_hi,int *x_r);
void gs_density(double *den, int *x_n_anal, int *n_theta, double *x_theta,
               double *I, double *a, double *b, double *x_z, int *z_len, int *x_r);
void std_nor_pts(int *r,double *bounds,double *z,double *w);

#ifdef __cplusplus
}
#endif

void h_1(double theta,int m,double *wgt,double I,double *z,double *h);
void h_update(double theta, double *wgt,
             int m1, double i_km1,double *zkm1,double *hkm1,
             int m2, double i_k,  double *zk,  double *hk);
void grid_pts_1(int r,double mu,double *z);
int grid_pts(int r,double mu,double a, double b, double *z,double *w);
int grid_pts_2(int r,double mu,double a, double b, double *z,double *w);
double prob_neg(double theta,int m,double ak,double *z,double *h,
               double i_km1,double i_k);
double prob_neg_2(double theta,int m,double ak,double *z,double *h,
                double i_km1,double i_k);
double prob_pos(double theta,int m,double bk,double *z,double *h,
               double i_km1,double i_k);
double prob_pos_2(double theta,int m,double bk,double *z,double *h,
                double i_km1,double i_k);

#endif
