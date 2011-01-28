#include "math.h"


int min (int a, int b)
{
  if (a < b)
    {
      return a;
    }
  else
    {
      return b;
    }
}


int max (int a, int b)
{
  if (a > b)
    {
      return a;
    }
  else
    {
      return b;
    }
}


int nint (double a)
{
  if (a > (int)(a + 0.5F))
    {
      return 1 + (int)a;
    }
  else
    {
      return (int)a;
    }
}


void Rotate (double x1, double y1, double z1,
	     double sn1, double cs1, double sn2, double cs2,
	     double *x2, double *y2, double *z2)
{
  double temp = z1 * cs2 - y1 * sn2;
  
  *x2 = temp * sn1 + x1 * cs1;
  *y2 =   z1 * sn2 + y1 * cs2;
  *z2 = temp * cs1 - x1 * sn1;
}


void AntiRotate (double x1, double y1, double z1,
		 double sn1, double cs1, double sn2, double cs2,
		 double *x2, double *y2, double *z2)
{
  double temp = cs1 * z1 + sn1 * x1;
  
  *x2 = cs1 * x1   - sn1 * z1;
  *y2 = cs2 * y1   - sn2 * temp;
  *z2 = cs2 * temp + sn2 * y1;
}


double ScalarProd (double x1[3], double x2[3])
{
  return x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
}


void VectorProd (double x1[3], double x2[3], double x3[3])
{
  x3[0] = x2[1] * x1[2] - x2[2] * x1[1];
  x3[1] = x2[2] * x1[0] - x2[0] * x1[2];
  x3[2] = x2[0] * x1[1] - x2[1] * x1[0];
}
