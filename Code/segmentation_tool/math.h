#ifndef MATH
#define MATH

#include <math.h>


int min (int a, int b);
int max (int a, int b);
int nint (double a);
void Rotate (double x1, double y1, double z1,
	     double sn1, double cs1, double sn2, double cs2,
	     double *x2, double *y2, double *z2);
void AntiRotate (double x1, double y1, double z1,
		 double sn1, double cs1, double sn2, double cs2,
		 double *x2, double *y2, double *z2);
double ScalarProd (double x1[3], double x2[3]);
void VectorProd (double x1[3], double x2[3], double x3[3]);


#endif // MATH
