#ifndef COLLISIONS_H_ //COLLISIONS_H_
#define COLLISIONS_H_ //COLLISIONS_H_
#include "net.h"

class Collision
{
  public:
    virtual void DoCollisions(double omega, int i, double *density, double *v_x,
                              double *v_y, double *v_z, double f_neq[], Net* net) = 0;

  protected:
    void DensityAndVelocity(double f[], double *density, double *v_x, double *v_y,
                            double *v_z);
    void CalculateFeq(double density, double v_x, double v_y, double v_z, double f_eq[]);
};

class MidFluidCollision : public Collision
{
};

class WallCollision : public Collision
{
};

class InletOutletCollision : public Collision
{
};

class InletOutletWallCollision : public Collision
{
};

class SimpleCollideAndStream : public MidFluidCollision
{
    void DoCollisions(double omega, int i, double *density, double *v_x, double *v_y,
                      double *v_z, double f_neq[], Net* net);
};

class ZeroVelocityEquilibrium : public WallCollision
{
    void DoCollisions(double omega, int i, double *density, double *v_x, double *v_y,
                      double *v_z, double f_neq[], Net* net);
};

class NonZeroVelocityBoundaryDensity : public InletOutletCollision
{
  public:
    NonZeroVelocityBoundaryDensity(double* iBounaryDensityArray);

    void DoCollisions(double omega, int i, double *density, double *v_x, double *v_y,
                      double *v_z, double f_neq[], Net* net);

  private:
    double* mBoundaryDensityArray;
};

class ZeroVelocitySetBoundaryDensity : public InletOutletWallCollision
{
  public:
    ZeroVelocitySetBoundaryDensity(double* iBoundaryDensityArray);
    void DoCollisions(double omega, int i, double *density, double *v_x, double *v_y,
                      double *v_z, double f_neq[], Net* net);

  private:
    double* mBoundaryDensityArray;
};

#endif //COLLISIONS_H_
