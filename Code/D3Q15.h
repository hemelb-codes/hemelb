#ifndef HEMELB_D3Q15_H
#define HEMELB_D3Q15_H

class D3Q15
{
  public:
    // The number of discrete velocity vectors
    static const int NUMVECTORS = 15;

    // The x, y and z components of each of the discrete velocity vectors
    static const int CX[NUMVECTORS];
    static const int CY[NUMVECTORS];
    static const int CZ[NUMVECTORS];

    // The index of the inverse direction of each discrete velocity vector
    static const int INVERSEDIRECTIONS[NUMVECTORS];

    // Functions to calculate the density and velocity from an array of fs and to calculate the
    // equilibrium f array from density and velocity.
    static void CalculateDensityAndVelocity(double f[], double *density, double *v_x,
      double *v_y, double *v_z);
    static void CalculateFeq(double density, double v_x, double v_y, double v_z,
      double f_eq[]);
    static void CalculateDensityVelocityFEq(double f[], double *density, double *v_x,
      double *v_y, double *v_z, double f_eq[]);
    static void CalculateVonMisesStress(double f[], double *stress,
      double iStressParameter);
    static void CalculateShearStress(double density, double f[], double nor[],
      double *stress, double iStressParameter);
};

#endif /* HEMELB_D3Q15_H */
