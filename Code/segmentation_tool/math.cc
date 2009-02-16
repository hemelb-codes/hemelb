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


int nint (float a)
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


void Rotate (float x1[3], float longitude, float latitude, float x2[3])
{
  float sn1, cs1;
  float sn2, cs2;
  float temp;
  
  
  sn1 = sinf(longitude);
  cs1 = cosf(longitude);
  sn2 = sinf(latitude);
  cs2 = cosf(latitude);

  temp = cs2 * x1[2] - sn2 * x1[1];

  x2[0] = sn1 * temp  + cs1 * x1[0];
  x2[1] = sn2 * x1[2] + cs2 * x1[1];
  x2[2] = cs1 * temp  - sn1 * x1[0];
}


void Rotate (float x1[3], float sn1, float cs1, float sn2, float cs2, float x2[3])
{
  float temp;
  
  
  temp = cs2 * x1[2] - sn2 * x1[1];

  x2[0] = sn1 * temp  + cs1 * x1[0];
  x2[1] = sn2 * x1[2] + cs2 * x1[1];
  x2[2] = cs1 * temp  - sn1 * x1[0];
}


void AntiRotate (float x1[3], float longitude, float latitude, float x2[3])
{
  float sn1, cs1;
  float sn2, cs2;
  float temp;
  
  
  sn1 = sinf(longitude);
  cs1 = cosf(longitude);
  sn2 = sinf(latitude);
  cs2 = cosf(latitude);
  
  temp = cs1 * x1[2] + sn1 * x1[0];

  x2[0] = cs1 * x1[0] - sn1 * x1[2];
  x2[1] = cs2 * x1[1] - sn2 * temp;
  x2[2] = cs2 * temp  + sn2 * x1[1];
}


void AntiRotate (float x1[3], float sn1, float cs1, float sn2, float cs2, float x2[3])
{
  float temp;
  
  
  temp = cs1 * x1[2] + sn1 * x1[0];
  
  x2[0] = cs1 * x1[0] - sn1 * x1[2];
  x2[1] = cs2 * x1[1] - sn2 * temp;
  x2[2] = cs2 * temp  + sn2 * x1[1];
}
