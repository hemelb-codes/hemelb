#ifndef HEMELB_D3Q15_H
#define HEMELB_D3Q15_H

#include "constants.h"

class D3Q15
{
  public:
    // The number of discrete velocity vectors
    static const unsigned int NUMVECTORS = 15;

    // The x, y and z components of each of the discrete velocity vectors
    static const int CX[NUMVECTORS];
    static const int CY[NUMVECTORS];
    static const int CZ[NUMVECTORS];

    // The index of the inverse direction of each discrete velocity vector
    static const int INVERSEDIRECTIONS[NUMVECTORS];

    // Functions to calculate the density and velocity from an array of SharedFCount and to calculate the
    // equilibrium f array from density and velocity.
    static void CalculateDensityAndVelocity(const distribn_t f[],
                                            distribn_t &density,
                                            distribn_t &v_x,
                                            distribn_t &v_y,
                                            distribn_t &v_z);
    static void CalculateFeq(const distribn_t &density,
                             const distribn_t &v_x,
                             const distribn_t &v_y,
                             const distribn_t &v_z,
                             distribn_t f_eq[]);
    static void CalculateDensityVelocityFEq(const distribn_t f[],
                                            distribn_t &density,
                                            distribn_t &v_x,
                                            distribn_t &v_y,
                                            distribn_t &v_z,
                                            distribn_t f_eq[]);
    static void CalculateVonMisesStress(const distribn_t f[],
                                        distribn_t &stress,
                                        const double iStressParameter);
    static void CalculateShearStress(const distribn_t &density,
                                     const distribn_t f[],
                                     const double nor[],
                                     distribn_t &stress,
                                     const double &iStressParameter);
};

#endif /* HEMELB_D3Q15_H */
