#include "gs_design_crt.h"

void std_nor_pts(int *r,double *bounds,double *z,double *w)
{  int grid_pts(int,double,double,double,double *,double *);
	grid_pts(r[0],0.,bounds[0],bounds[1],z,w);
}
