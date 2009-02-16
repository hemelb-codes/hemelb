#ifndef MATH
#define MATH

#include <math.h>


int min (int a, int b);
int max (int a, int b);
int nint (float a);
void Rotate (float x1[3], float longitude, float latitude, float x2[3]);
void Rotate (float x1[3], float sn1, float cs1, float sn2, float cs2, float x2[3]);
void AntiRotate (float x1[3], float longitude, float latitude, float x2[3]);
void AntiRotate (float x1[3], float sn1, float cs1, float sn2, float cs2, float x2[3]);


#endif // MATH
