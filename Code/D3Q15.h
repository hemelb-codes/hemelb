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
    static void CalculateDensityAndVelocity(const double f[],
                                            double &density,
                                            double &v_x,
                                            double &v_y,
                                            double &v_z);
    static void CalculateFeq(const double &density,
                             const double &v_x,
                             const double &v_y,
                             const double &v_z,
                             double f_eq[]);
    static void CalculateDensityVelocityFEq(const double f[],
                                            double &density,
                                            double &v_x,
                                            double &v_y,
                                            double &v_z,
                                            double f_eq[]);
    static void CalculateVonMisesStress(const double f[],
                                        double &stress,
                                        const double iStressParameter);
    static void CalculateShearStress(const double &density,
                                     const double f[],
                                     const double nor[],
                                     double &stress,
                                     const double &iStressParameter);
};

#endif /* HEMELB_D3Q15_H */
