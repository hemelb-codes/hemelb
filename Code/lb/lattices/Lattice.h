#ifndef HEMELB_LB_LATTICES_LATTICE_H
#define HEMELB_LB_LATTICES_LATTICE_H

#include <cmath>
#include "constants.h"
#include "lb/lattices/LatticeInfo.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"

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

  namespace lb
  {
    namespace lattices
    {
      template<class DmQn>
      class Lattice
      {
        public:
          static void CalculateDensityAndVelocity(const distribn_t f[],
                                                  distribn_t &density,
                                                  distribn_t &v_x,
                                                  distribn_t &v_y,
                                                  distribn_t &v_z)
          {
            density = v_x = v_y = v_z = 0.0;

            for (Direction direction = 0; direction < DmQn::NUMVECTORS; ++direction)
            {
              density += f[direction];
              v_x += DmQn::CX[direction] * f[direction];
              v_y += DmQn::CY[direction] * f[direction];
              v_z += DmQn::CZ[direction] * f[direction];
            }
          }

          static void CalculateFeq(const distribn_t &density,
                                   const distribn_t &v_x,
                                   const distribn_t &v_y,
                                   const distribn_t &v_z,
                                   distribn_t f_eq[])
          {
            const distribn_t density_1 = 1. / density;
            const distribn_t velocityMagnitudeSquared = v_x * v_x + v_y * v_y + v_z * v_z;

            for (Direction direction = 0; direction < DmQn::NUMVECTORS; ++direction)
            {
              const distribn_t velocityComponentInThisDirection = DmQn::CX[direction] * v_x
                  + DmQn::CY[direction] * v_y + DmQn::CZ[direction] * v_z;

              f_eq[direction] = DmQn::EQMWEIGHTS[direction] * (density - (3. / 2.)
                  * velocityMagnitudeSquared * density_1 + (9. / 2.) * density_1
                  * velocityComponentInThisDirection * velocityComponentInThisDirection + 3.
                  * velocityComponentInThisDirection);
            }
          }

          // Calculate density, velocity and the equilibrium distribution
          // functions according to the D3Q15 model.  The calculated v_x, v_y
          // and v_z are actually density * velocity, because we are using the
          // compressible model.
          static void CalculateDensityVelocityFEq(const distribn_t f[],
                                                  distribn_t &density,
                                                  distribn_t &v_x,
                                                  distribn_t &v_y,
                                                  distribn_t &v_z,
                                                  distribn_t f_eq[])
          {
            CalculateDensityAndVelocity(f, density, v_x, v_y, v_z);

            CalculateFeq(density, v_x, v_y, v_z, f_eq);
          }

          static void CalculateEntropicDensityVelocityFEq(const distribn_t f[],
                                                          distribn_t &density,
                                                          distribn_t &v_x,
                                                          distribn_t &v_y,
                                                          distribn_t &v_z,
                                                          distribn_t f_eq[])
          {
            CalculateDensityAndVelocity(f, density, v_x, v_y, v_z);

            CalculateEntropicFeq(density, v_x, v_y, v_z, f_eq);
          }

          // von Mises stress computation given the non-equilibrium distribution functions.
          static void CalculateVonMisesStress(const distribn_t f[],
                                              distribn_t &stress,
                                              const double iStressParameter)
          {
            // Recall that sigma_ij = Sum_l f(l) * C_il * C_jl

            // First calculate sigma_xx - sigma_yy.
            // Using standard form of sigma_ij, sigma_xx - sigma_yy
            //   = Sum_l f(l) * (Cx(l) * Cx(l) - Cy(l) * Cy(l))
            // We calculate sigma_yy - sigma_zz and sigma_xx and sigma_zz
            // in the same way.
            distribn_t sigma_xx_yy = 0.0;
            distribn_t sigma_yy_zz = 0.0;
            distribn_t sigma_xx_zz = 0.0;

            // We will also require sigma_xy, sigma_yz and sigma_xz, calculated in the usual
            // way.
            distribn_t sigma_xy = 0.0;
            distribn_t sigma_xz = 0.0;
            distribn_t sigma_yz = 0.0;

            for (Direction direction = 0; direction < DmQn::NUMVECTORS; ++direction)
            {
              sigma_xx_yy += f[direction] * (DmQn::CX[direction] * DmQn::CX[direction]
                  - DmQn::CY[direction] * DmQn::CY[direction]);
              sigma_yy_zz += f[direction] * (DmQn::CY[direction] * DmQn::CY[direction]
                  - DmQn::CZ[direction] * DmQn::CZ[direction]);
              sigma_xx_zz += f[direction] * (DmQn::CX[direction] * DmQn::CX[direction]
                  - DmQn::CZ[direction] * DmQn::CZ[direction]);

              sigma_xy += f[direction] * DmQn::CX[direction] * DmQn::CY[direction];
              sigma_xz += f[direction] * DmQn::CX[direction] * DmQn::CZ[direction];
              sigma_yz += f[direction] * DmQn::CY[direction] * DmQn::CZ[direction];
            }

            distribn_t a = sigma_xx_yy * sigma_xx_yy + sigma_yy_zz * sigma_yy_zz + sigma_xx_zz
                * sigma_xx_zz;
            distribn_t b = sigma_xy * sigma_xy + sigma_xz * sigma_xz + sigma_yz * sigma_yz;

            stress = iStressParameter * sqrt(a + 6.0 * b);
          }

          // The magnitude of the tangential component of the shear stress acting on the
          // wall.
          static void CalculateShearStress(const distribn_t &density,
                                           const distribn_t f[],
                                           const util::Vector3D<double> nor,
                                           distribn_t &stress,
                                           const double &iStressParameter)
          {
            // sigma_ij is the force
            // per unit area in
            // direction i on the
            // plane with the normal
            // in direction j
            distribn_t stress_vector[] = { 0.0, 0.0, 0.0 }; // Force per unit area in
            // direction i on the
            // plane perpendicular to
            // the surface normal
            distribn_t square_stress_vector = 0.0;
            distribn_t normal_stress = 0.0; // Magnitude of force per
            // unit area normal to the
            // surface

            distribn_t temp = iStressParameter * (-sqrt(2.0));

            Order2Tensor pi = CalculatePiTensor(f);

            for (unsigned i = 0; i < 3; i++)
            {
              for (unsigned j = 0; j < 3; j++)
                stress_vector[i] += pi[i][j] * nor[j] * temp;

              square_stress_vector += stress_vector[i] * stress_vector[i];
              normal_stress += stress_vector[i] * nor[i];
            }
            // shear_stress^2 + normal_stress^2 = stress_vector^2
            stress = sqrt(square_stress_vector - normal_stress * normal_stress);
          }

          static Order2Tensor CalculatePiTensor(const distribn_t* const f)
          {
            Order2Tensor ret;

            // Fill in 0,0 1,0 1,1 2,0 2,1 2,2
            for (int ii = 0; ii < 3; ++ii)
            {
              for (int jj = 0; jj <= ii; ++jj)
              {
                ret[ii][jj] = 0.0;
                for (unsigned int l = 0; l < DmQn::NUMVECTORS; ++l)
                {
                  ret[ii][jj] += f[l] * DmQn::discreteVelocityVectors[ii][l]
                      * DmQn::discreteVelocityVectors[jj][l];
                }
              }
            }

            // Exploit the symmetry to fill in 0,1 0,2 1,2
            for (int ii = 0; ii < 3; ++ii)
            {
              for (int jj = ii + 1; jj < 3; ++jj)
              {
                ret[ii][jj] = ret[jj][ii];
              }
            }

            return ret;
          }

          static distribn_t CalculateShearRate(const distribn_t &iTau,
                                               const distribn_t iFNeq[],
                                               const distribn_t &iDensity)
          {
            distribn_t shear_rate = 0.0;
            distribn_t strain_rate_tensor_i_j;

            for (unsigned row = 0; row < 3; row++)
            {
              for (unsigned column = 0; column < 3; column++)
              {
                strain_rate_tensor_i_j = CalculateStrainRateTensorComponent(row,
                                                                            column,
                                                                            iTau,
                                                                            iFNeq,
                                                                            iDensity);
                shear_rate += strain_rate_tensor_i_j * strain_rate_tensor_i_j;
              }
            }

            shear_rate = sqrt(shear_rate);

            return shear_rate;
          }

          // Entropic ELBM has an analytical form for FEq
          // (see Aidun and Clausen "Lattice-Boltzmann Method for Complex Flows" Annu. Rev. Fluid. Mech. 2010)
          static void CalculateEntropicFeq(const distribn_t &density,
                                           const distribn_t &v_x,
                                           const distribn_t &v_y,
                                           const distribn_t &v_z,
                                           distribn_t f_eq[])
          {
            // Get velocity
            util::Vector3D<distribn_t> velocity = util::Vector3D<distribn_t>(v_x, v_y, v_z)
                / density;

            // Combining some terms for use in evaluating the next few terms
            // B_i = sqrt(1 + 3 * u_i^2)
            util::Vector3D<distribn_t> B = util::Vector3D<distribn_t>(sqrt(1.0 + 3.0 * velocity.x
                                                                          * velocity.x),
                                                                      sqrt(1.0 + 3.0 * velocity.y
                                                                          * velocity.y),
                                                                      sqrt(1.0 + 3.0 * velocity.z
                                                                          * velocity.z));

            // The formula has contains the product term1_i*(term2_i)^e_ia
            // term1_i is 2 - B_i
            util::Vector3D<distribn_t> term1 = util::Vector3D<distribn_t>(2.0) - B;

            // term2_i is (2*u_i + B)/(1 - u_i)
            util::Vector3D<distribn_t> term2 =
                (velocity * 2.0 + B).PointwiseDivision(util::Vector3D<distribn_t>::Unity()
                    - velocity);

            for (Direction direction = 0; direction < DmQn::NUMVECTORS; ++direction)
            {
              f_eq[direction] = density * DmQn::EQMWEIGHTS[direction] * term1.x * term1.y * term1.z
                  * util::NumericalFunctions::IntegerPower(term2.x, DmQn::CX[direction])
                  * util::NumericalFunctions::IntegerPower(term2.y, DmQn::CY[direction])
                  * util::NumericalFunctions::IntegerPower(term2.z, DmQn::CZ[direction]);
            }
          }

          static LatticeInfo* GetLatticeInfo()
          {
            if (singletonInfo == NULL)
            {
              util::Vector3D<int> vectors[DmQn::NUMVECTORS];
              int inverseVectorIndices[DmQn::NUMVECTORS];

              for (Direction direction = 0; direction < DmQn::NUMVECTORS; ++direction)
              {
                vectors[direction] = util::Vector3D<int>(DmQn::CX[direction],
                                                         DmQn::CY[direction],
                                                         DmQn::CZ[direction]);
                inverseVectorIndices[direction] = DmQn::INVERSEDIRECTIONS[direction];
              }

              singletonInfo = new LatticeInfo(DmQn::NUMVECTORS, vectors, inverseVectorIndices);
            }

            return singletonInfo;
          }

        private:
          static distribn_t CalculateStrainRateTensorComponent(const unsigned &iRow,
                                                               const unsigned &iColumn,
                                                               const distribn_t &iTau,
                                                               const distribn_t iFNeq[],
                                                               const distribn_t &iDensity)
          {
            distribn_t strain_rate_tensor_i_j = 0.0;

            for (unsigned vec_index = 0; vec_index < DmQn::NUMVECTORS; vec_index++)
            {
              strain_rate_tensor_i_j += iFNeq[vec_index]
                  * (DmQn::discreteVelocityVectors[iRow][vec_index]
                      * DmQn::discreteVelocityVectors[iColumn][vec_index]);
            }

            strain_rate_tensor_i_j *= -1.0 / (2.0 * iTau * iDensity * Cs2);

            return strain_rate_tensor_i_j;
          }

          static LatticeInfo* singletonInfo;
      };
    }
  }
}

#endif
