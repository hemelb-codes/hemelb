#ifndef HEMELB_D3Q15_H
#define HEMELB_D3Q15_H

#include "constants.h"
namespace hemelb
{
  struct Order2Tensor
  {
      distribn_t tensor[3][3];

      /**
       * Convenience accessor.
       *
       * * @param row
       * @return
       */
      distribn_t* operator [](const unsigned int row)
      {
        return tensor[row];
      }
  };

  class D3Q15
  {
    public:
      // The number of discrete velocity vectors
      static const unsigned int NUMVECTORS = 15;

      // The x, y and z components of each of the discrete velocity vectors
      static const int CX[NUMVECTORS];
      static const int CY[NUMVECTORS];
      static const int CZ[NUMVECTORS];
      static const int* discreteVelocityVectors[3];

      static const double EQMWEIGHTS[NUMVECTORS];

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
      static void CalculateEntropicFeq(const distribn_t &density,
                                       const distribn_t &v_x,
                                       const distribn_t &v_y,
                                       const distribn_t &v_z,
                                       distribn_t f_eq[]);
      static void CalculateEntropicDensityVelocityFEq(const distribn_t f[],
                                                      distribn_t &density,
                                                      distribn_t &v_x,
                                                      distribn_t &v_y,
                                                      distribn_t &v_z,
                                                      distribn_t f_eq[]);
      static void CalculateVonMisesStress(const distribn_t f[],
                                          distribn_t &stress,
                                          const double iStressParameter);

      /**
       * Computes the Pi tensor, the second-order moment of the f-distribution.
       * Pi = Sum over directions {c_i c_i f_i}
       *
       * @param f
       * @return
       */
      static Order2Tensor CalculatePiTensor(const distribn_t* const f);

      static void CalculateShearStress(const distribn_t &density,
                                       const distribn_t f[],
                                       const double nor[],
                                       distribn_t &stress,
                                       const double &iStressParameter);
      /*
       * Computes the (j,k)-th entry of the strain rate tensor according to the expression
       *
       *    S = (-1 / (2 * tau * rho * Cs2)) * sum_over_i (e_i*e_i*f^(neq)[i])
       *
       * where e_i is the i-th direction in the lattice, and e_i*e_i is an outter product.
       *
       * @param iRow tensor row.
       * @param iColumn tensor column.
       * @param iTau relaxation time (dimensionless).
       * @param iFNeq non-equilibrium velocity distribution function.
       * @param iDensity local density (dimensionless).
       *
       * @return S[iRow,iColumn] (dimensionless).
       */
      static distribn_t CalculateStrainRateTensorComponent(const unsigned &iRow,
                                                           const unsigned &iColumn,
                                                           const distribn_t &iTau,
                                                           const distribn_t iFNeq[],
                                                           const distribn_t &iDensity);
      /*
       * Computes shear-rate according to the expression
       *
       *    gamma_dot = sqrt( sum_over_i_j( S[i,j]*S[i,j] ) )
       *
       * where S is the strain-rate tensor.
       *
       * TODO Different expressions for the shear-rate can be found in different
       * publications: 2*sqrt(s_ij * s_ij), sqrt(2 * s_ij * s_ij), sqrt(s_ij * s_ij)
       *
       * @param iTau relaxation time (dimensionless).
       * @param iFNeq non-equilibrium velocity distribution function.
       * @param iDensity local density (dimensionless).
       *
       * @return gamma (dimensionless)
       */
      static distribn_t CalculateShearRate(const distribn_t &iTau,
                                           const distribn_t iFNeq[],
                                           const distribn_t &iDensity);
  };
}
#endif /* HEMELB_D3Q15_H */
