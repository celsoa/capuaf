#include<math.h>
#include<stdio.h>
#include<stdlib.h>  // drand48
#include<float.h> 
#include "cap.h"

void u2beta_vec(float *pArray_u, float *pArray_beta, int nsol);
float beta2delta(float beta);
void beta2delta_vec(float *pArray_beta, float *pArray_delta, int nsol);
float gamma2v(float gamma);
float v2gamma(float v);
void v2gamma_vec(float *pArray_v, int nsol, float *pArray_gamma);
float dip2h(float dip);
float h2dip(float h);
void h2dip_vec(float *pArray_h, int nsol, float *pArray_dip);
float siso2delta(float siso);
void siso2delta_vec(float *pArray_s, int nsol, float *pArray_dip);
void randvec(float xmin, float xmax, int nsol, float *pArray);
void gridvec(float xmin, float xmax, int npoints, float *pArray);
void magvec(float xmin, float xmax, float dx, float *pArray);
float delta2beta(float delta);
float w2u(float w);
float u2w(float u);
float delta2w(float delta);
int find_nearest_index(float value, float * pData, int len);
void interp_lin(float *x, float *y, int x_size, float *xx, float *yy, int xx_size);
