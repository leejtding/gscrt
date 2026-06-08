#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "gs_design_crt.h"

/* Register all .C entry points in gscrt. */

/*
void gs_bound(int *x_n_anal,double *I,double *a,double *b,double *prob_lo,double *prob_hi,
             double *x_tol,int *x_r,int *retval,int *printerr) */

static R_NativePrimitiveArgType gsbound_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gs_upper_1(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
              double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr) */

static R_NativePrimitiveArgType gs_upper_1_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gs_upper_2(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
             double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr) */

static R_NativePrimitiveArgType gs_upper_2_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gs_lower_1(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
              double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr) */

static R_NativePrimitiveArgType gs_lower_1_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gs_lower_2(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
              double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr) */

static R_NativePrimitiveArgType gs_lower_2_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gs_bounds_1(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
               double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr) */

static R_NativePrimitiveArgType gs_bounds_1_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gs_bounds_2(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
               double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr) */

static R_NativePrimitiveArgType gs_bounds_2_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gs_bounds_nb_1(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
               double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr) */

static R_NativePrimitiveArgType gs_bounds_nb_1_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void gs_bounds_nb_2(int *x_n_anal,double *x_theta,double *I,double *a,double *b,double *prob_lo,
               double *prob_hi,double *x_tol,int *x_r,int *retval,int *printerr) */

static R_NativePrimitiveArgType gs_bounds_nb_2_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/*
void prob_rej_1(int *x_n_anal,int *n_theta,double *x_theta,double *I,double *a,double *b,
              double *x_prob_lo,double *x_prob_hi,int *x_r) */

static R_NativePrimitiveArgType prob_rej_1_t[] = {
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

/*
void prob_rej_2(int *x_n_anal,int *n_theta,double *x_theta,double *I,double *a,double *b,
              double *x_prob_lo,double *x_prob_hi,int *x_r) */

static R_NativePrimitiveArgType prob_rej_2_t[] = {
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

/*
void prob_stop_1(int *x_n_anal,int *n_theta,double *x_theta,double *I,double *a,double *b,
               double *x_prob_lo,double *x_prob_hi,int *x_r) */

static R_NativePrimitiveArgType prob_stop_1_t[] = {
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

/*
void prob_stop_2(int *x_n_anal,int *n_theta,double *x_theta,double *I,double *a,double *b,
               double *x_prob_lo,double *x_prob_hi,int *x_r) */

static R_NativePrimitiveArgType prob_stop_2_t[] = {
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

/*
void gs_density(double *den, int *x_n_anal, int *n_theta, double *x_theta,
               double *I, double *a, double *b, double *x_z,
               int *z_len, int *x_r) */

static R_NativePrimitiveArgType gs_density_t[] = {
  REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP
};

/*
void std_nor_pts(int *r,double *bounds,double *z,double *w) */

static R_NativePrimitiveArgType std_nor_pts_t[] = {
  INTSXP, REALSXP, REALSXP, REALSXP
};

/* define array of all C entry points */
static const R_CMethodDef c_entries[] = {
  {"gs_bound", (DL_FUNC) &gs_bound, 10, gsbound_t},
  {"gs_upper_1", (DL_FUNC) &gs_upper_1, 11, gs_upper_1_t},
  {"gs_upper_2", (DL_FUNC) &gs_upper_2, 11, gs_upper_2_t},
  {"gs_lower_1", (DL_FUNC) &gs_lower_1, 11, gs_lower_1_t},
  {"gs_lower_2", (DL_FUNC) &gs_lower_2, 11, gs_lower_2_t},
  {"gs_bounds_1", (DL_FUNC) &gs_bounds_1, 11, gs_bounds_1_t},
  {"gs_bounds_2", (DL_FUNC) &gs_bounds_2, 11, gs_bounds_2_t},
  {"gs_bounds_nb_1", (DL_FUNC) &gs_bounds_nb_1, 11, gs_bounds_nb_1_t},
  {"gs_bounds_nb_2", (DL_FUNC) &gs_bounds_nb_2, 11, gs_bounds_nb_2_t},
  {"prob_rej_1", (DL_FUNC) &prob_rej_1, 9, prob_rej_1_t},
  {"prob_rej_2", (DL_FUNC) &prob_rej_2, 9, prob_rej_2_t},
  {"prob_stop_1", (DL_FUNC) &prob_stop_1, 9, prob_stop_1_t},
  {"prob_stop_2", (DL_FUNC) &prob_stop_2, 9, prob_stop_2_t},
  {"gs_density", (DL_FUNC) &gs_density, 10, gs_density_t},
  {"std_nor_pts", (DL_FUNC) &std_nor_pts, 4, std_nor_pts_t},
  {NULL, NULL, 0, NULL}
};

/* now register in the init function */
extern "C" void R_init_gscrt(DllInfo *dll)
{
  R_registerRoutines(dll, c_entries, NULL, NULL, NULL);

  /* Resolve .C calls through the registration table. */
  R_useDynamicSymbols(dll, FALSE);

  /* allow .C calls by character strings: */
  R_forceSymbols(dll, FALSE);
}
