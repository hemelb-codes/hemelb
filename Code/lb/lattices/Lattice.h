// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_LATTICES_LATTICE_H
#define HEMELB_LB_LATTICES_LATTICE_H

#include <cmath>
#include "constants.h"
#include "lb/lattices/LatticeInfo.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"
#include "util/Matrix3D.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      template<class DmQn>
      class Lattice
      {
        public:
          inline static void CalculateDensityAndMomentum(const distribn_t f[],
                                                         distribn_t &density,
                                                         distribn_t &momentum_x,
                                                         distribn_t &momentum_y,
                                                         distribn_t &momentum_z)
          {
            density = momentum_x = momentum_y = momentum_z = 0.0;

            for (Direction direction = 0; direction < DmQn::NUMVECTORS; ++direction)
            {
              density += f[direction];
              momentum_x += DmQn::CX[direction] * f[direction];
              momentum_y += DmQn::CY[direction] * f[direction];
              momentum_z += DmQn::CZ[direction] * f[direction];
            }
          }

          inline static void CalculateFeq(const distribn_t &density,
                                          const distribn_t &momentum_x,
                                          const distribn_t &momentum_y,
                                          const distribn_t &momentum_z,
                                          distribn_t f_eq[])
          {
            const distribn_t density_1 = 1. / density;
            const distribn_t momentumMagnitudeSquared = momentum_x * momentum_x + momentum_y * momentum_y
                + momentum_z * momentum_z;

            for (Direction direction = 0; direction < DmQn::NUMVECTORS; ++direction)
            {
              const distribn_t momentumComponentInThisDirection = DmQn::CX[direction] * momentum_x
                  + DmQn::CY[direction] * momentum_y + DmQn::CZ[direction] * momentum_z;

              f_eq[direction] = DmQn::EQMWEIGHTS[direction]
                  * (density - (3. / 2.) * momentumMagnitudeSquared * density_1
                      + (9. / 2.) * density_1 * momentumComponentInThisDirection * momentumComponentInThisDirection
                      + 3. * momentumComponentInThisDirection);
            }
          }

          // Calculate density, momentum and the equilibrium distribution
          // functions according to the D3Q15 model.  The calculated momentum_x, momentum_y
          // and momentum_z are actually density * velocity, because we are using the
          // compressible model.
          inline static void CalculateDensityMomentumFEq(const distribn_t f[],
                                                         distribn_t &density,
                                                         distribn_t &momentum_x,
                                                         distribn_t &momentum_y,
                                                         distribn_t &momentum_z,
                                                         distribn_t f_eq[])
          {
            CalculateDensityAndMomentum(f, density, momentum_x, momentum_y, momentum_z);

            CalculateFeq(density, momentum_x, momentum_y, momentum_z, f_eq);
          }

          // von Mises stress computation given the non-equilibrium distribution functions.
          inline static void CalculateVonMisesStress(const distribn_t f[],
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
              sigma_xx_yy += f[direction]
                  * (DmQn::CX[direction] * DmQn::CX[direction] - DmQn::CY[direction] * DmQn::CY[direction]);
              sigma_yy_zz += f[direction]
                  * (DmQn::CY[direction] * DmQn::CY[direction] - DmQn::CZ[direction] * DmQn::CZ[direction]);
              sigma_xx_zz += f[direction]
                  * (DmQn::CX[direction] * DmQn::CX[direction] - DmQn::CZ[direction] * DmQn::CZ[direction]);

              sigma_xy += f[direction] * DmQn::CX[direction] * DmQn::CY[direction];
              sigma_xz += f[direction] * DmQn::CX[direction] * DmQn::CZ[direction];
              sigma_yz += f[direction] * DmQn::CY[direction] * DmQn::CZ[direction];
            }

            distribn_t a = sigma_xx_yy * sigma_xx_yy + sigma_yy_zz * sigma_yy_zz + sigma_xx_zz * sigma_xx_zz;
            distribn_t b = sigma_xy * sigma_xy + sigma_xz * sigma_xz + sigma_yz * sigma_yz;

            stress = iStressParameter * sqrt(a + 6.0 * b);
          }

          /**
           * Calculates the traction vector on a surface point (units of stress). This is done by multiplying the full
           * stress tensor by the (outward pointing) surface normal at that point.
           *
           *    \vec{t} = \sigma \dot \vec{normal}
           *
           * @param density density at a given site
           * @param tau relaxation time
           * @param fNonEquilibrium non equilibrium part of the distribution function
           * @param wallNormal wall normal at a given point
           * @param tractionVector traction vector at a given point
           */
          inline static void CalculateTractionVectorOnAPoint(const distribn_t density,
                                                             const distribn_t tau,
                                                             const distribn_t fNonEquilibrium[],
                                                             const util::Vector3D<DimensionlessQuantity>& wallNormal,
                                                             util::Vector3D<LatticeStress>& tractionVector)
          {
            util::Matrix3D sigma;
            CalculateStressTensor(density, tau, fNonEquilibrium, sigma);

            // Multiply the stress tensor by the surface normal
            sigma.timesVector(wallNormal, tractionVector);
          }

          /**
           * Calculates the projection of the traction vector on the plane tangential to the geometry surface (defined by current surface
           * point and the normal vector provided). This is done with the following formula:
           *
           *    \vec{t_tan} = \vec{t} - dot(\vec{t}, \vec{n})*\vec{n}
           *
           * where t is the traction vector (see CalculateTractionVectorOnAPoint for definition) and n is the normal
           *
           * @param density density at a given site
           * @param tau relaxation time
           * @param fNonEquilibrium non equilibrium part of the distribution function
           * @param wallNormal wall normal at a given point
           * @param tractionTangentialComponent tangential projection of the traction vector
           */
          inline static void CalculateTangentialProjectionTractionVector(const distribn_t density,
                                                                         const distribn_t tau,
                                                                         const distribn_t fNonEquilibrium[],
                                                                         const util::Vector3D<DimensionlessQuantity>& wallNormal,
                                                                         util::Vector3D<LatticeStress>& tractionTangentialComponent)
          {
            util::Vector3D<LatticeStress> tractionVector;
            CalculateTractionVectorOnAPoint(density, tau, fNonEquilibrium, wallNormal, tractionVector);

            LatticeStress magnitudeNormalProjectionTraction = tractionVector.Dot(wallNormal);

            tractionTangentialComponent = tractionVector - wallNormal * magnitudeNormalProjectionTraction;
          }

          /**
           * Calculate the full stress tensor at a given fluid site (including both pressure and deviatoric part)
           *
           * The stress tensor is assembled based on the formula:
           *
           *    \sigma = p*I + 2*\mu*S = p*I - \Pi^{(neq)}
           *
           * where p is hydrostatic pressure, I is the identity tensor, S is the strain rate tensor, and \mu is the
           * viscosity. -2*\mu*S can be shown to be equals to the non equilibrium part of the moment flux tensor \Pi^{(neq)}.
           *
           * \Pi^{(neq)} is assumed to be defined as in Chen&Doolen 1998:
           *
           *    \Pi^{(neq)} = (1 - 1/(2*\tau)) * \sum_over_i e_i e_i f^{(neq)}_i
           *
           * where \tau is the relaxation time and e_i is the i-th direction vector
           *
           * @param density density at a given site
           * @param tau relaxation time
           * @param fNonEquilibrium non equilibrium part of the distribution function
           * @param stressTensor full stress tensor at a given site
           */
          inline static void CalculateStressTensor(const distribn_t density,
                                                   const distribn_t tau,
                                                   const distribn_t fNonEquilibrium[],
                                                   util::Matrix3D& stressTensor)
          {
            // Initialises the stress tensor to the deviatoric part, i.e. -\Pi^{(neq)}
            stressTensor = CalculatePiTensor(fNonEquilibrium);
            stressTensor *= 1 - 1 / (2 * tau);

            // Adds the pressure component to the stress tensor
            LatticePressure pressure = density * Cs2;
            stressTensor.addDiagonal(pressure);
          }

          /**
           * The magnitude of the tangential component of the shear stress acting on the
           * wall (i.e. tangential component of the traction vector). For this method to
           * make sense f has to be the non equilibrium part of a distribution function.
           *
           * The stress tensor computed in this method only includes the deviatoric part
           * and not the component corresponding to the pressure. Do not use the intermediate
           * traction value unless you understand the implications (use CalculateTractionVectorOnAPoint
           * instead).
           */
          inline static void CalculateWallShearStressMagnitude(const distribn_t &density,
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

            // Multiplying the second moment of the non equilibrium function by temp gives the non equilibrium part
            // of the moment flux tensor pi.
            distribn_t temp = iStressParameter * (-sqrt(2.0));

            // Computes the second moment of the equilibrium function f.
            util::Matrix3D pi = CalculatePiTensor(f);

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

          /**
           * Despite its name, this method does not compute the whole pi tensor (i.e. momentum flux tensor). What it does is
           * computing the second moment of a distribution function. If this distribution happens to be f_eq, the resulting
           * tensor will be the equilibrium part of pi. However, if the distribution function is f_neq, the result WON'T be
           * the non equilibrium part of pi. In order to get it, you will have to multiply by (1 - timestep/2*tau)
           *
           * @param f distribution function
           * @return second moment of the distribution function f
           */
          inline static util::Matrix3D CalculatePiTensor(const distribn_t* const f)
          {
            util::Matrix3D ret;

            // Fill in 0,0 1,0 1,1 2,0 2,1 2,2
            for (int ii = 0; ii < 3; ++ii)
            {
              for (int jj = 0; jj <= ii; ++jj)
              {
                ret[ii][jj] = 0.0;
                for (unsigned int l = 0; l < DmQn::NUMVECTORS; ++l)
                {
                  ret[ii][jj] += f[l] * DmQn::discreteVelocityVectors[ii][l] * DmQn::discreteVelocityVectors[jj][l];
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

          inline static distribn_t CalculateShearRate(const distribn_t &iTau,
                                                      const distribn_t iFNeq[],
                                                      const distribn_t &iDensity)
          {
            distribn_t shear_rate = 0.0;
            distribn_t strain_rate_tensor_i_j;

            for (unsigned row = 0; row < 3; row++)
            {
              for (unsigned column = 0; column < 3; column++)
              {
                strain_rate_tensor_i_j = CalculateStrainRateTensorComponent(row, column, iTau, iFNeq, iDensity);
                shear_rate += strain_rate_tensor_i_j * strain_rate_tensor_i_j;
              }
            }

            shear_rate = sqrt(shear_rate);

            return shear_rate;
          }

          // Entropic ELBM has an analytical form for FEq
          // (see Aidun and Clausen "Lattice-Boltzmann Method for Complex Flows" Annu. Rev. Fluid. Mech. 2010)
          // Originally Ansumali, S., Karlin, I. V., and Ottinger, H.C. (2003) Minimal entropic kinetic models
          // for hydrodynamics. Europhys. Lett. 63(6), 798–804
          inline static void CalculateEntropicFeqAnsumali(const distribn_t &density,
                                                          const distribn_t &momentum_x,
                                                          const distribn_t &momentum_y,
                                                          const distribn_t &momentum_z,
                                                          distribn_t f_eq[])
          {
            // Get velocity
            util::Vector3D<distribn_t> velocity = util::Vector3D<distribn_t>(momentum_x, momentum_y, momentum_z)
                / density;

            // Combining some terms for use in evaluating the next few terms
            // B_i = sqrt(1 + 3 * u_i^2)
            util::Vector3D<distribn_t> B = util::Vector3D<distribn_t>(sqrt(1.0 + 3.0 * velocity.x * velocity.x),
                                                                      sqrt(1.0 + 3.0 * velocity.y * velocity.y),
                                                                      sqrt(1.0 + 3.0 * velocity.z * velocity.z));

            // The formula contains the product term1_i*(term2_i)^e_ia
            // term1_i is 2 - B_i
            util::Vector3D<distribn_t> term1 = util::Vector3D<distribn_t>(2.0) - B;

            // term2_i is (2*u_i + B)/(1 - u_i)
            util::Vector3D<distribn_t> term2 = (velocity * 2.0 + B).PointwiseDivision(util::Vector3D<distribn_t>::Ones()
                - velocity);

            for (Direction direction = 0; direction < DmQn::NUMVECTORS; ++direction)
            {
              f_eq[direction] = density * DmQn::EQMWEIGHTS[direction] * term1.x * term1.y * term1.z
                  * util::NumericalFunctions::IntegerPower(term2.x, DmQn::CX[direction])
                  * util::NumericalFunctions::IntegerPower(term2.y, DmQn::CY[direction])
                  * util::NumericalFunctions::IntegerPower(term2.z, DmQn::CZ[direction]);
            }
          }

          /**
           * Calculate entropic equilibrium distribution, as in Chikatamarla et al (PRL, 97, 010201 (2006)
           *
           * @param density
           * @param momentum_x
           * @param momentum_y
           * @param momentum_z
           * @param f_eq
           */
          inline static void CalculateEntropicFeqChik(const distribn_t &density,
                                                      const distribn_t &momentum_x,
                                                      const distribn_t &momentum_y,
                                                      const distribn_t &momentum_z,
                                                      distribn_t f_eq[])
          {
            // Get velocity and the vector with velocity components squared.
            util::Vector3D<distribn_t> velocity = util::Vector3D<distribn_t>(momentum_x, momentum_y, momentum_z)
                / (density);
            util::Vector3D<distribn_t> velocitySquared = velocity.PointwiseMultiplication(velocity);
            util::Vector3D<distribn_t> velocityFour = velocitySquared.PointwiseMultiplication(velocitySquared);
            util::Vector3D<distribn_t> velocityEight = velocityFour.PointwiseMultiplication(velocityFour);

            // Compute in advance the first four powers of the velocity magnitude squared.
            distribn_t velocityMagnitudeSquared = velocity.GetMagnitudeSquared();
            distribn_t velocityMagnitudeFour = velocityMagnitudeSquared * velocityMagnitudeSquared;
            distribn_t velocityMagnitudeSix = velocityMagnitudeFour * velocityMagnitudeSquared;
            distribn_t velocityMagnitudeEight = velocityMagnitudeSix * velocityMagnitudeSquared;

            // Compute chi as per equation (9).
            distribn_t chi = 1.0 + (-3.0 * velocityMagnitudeSquared / 2.0) + 9.0 * velocityMagnitudeFour / 8.0;

            // Add in the (6) term.
            chi += 27.0
                * ( (-velocityMagnitudeSix)
                    + 2.0 * (velocitySquared.y + velocitySquared.z)
                        * (velocityMagnitudeSquared * velocitySquared.x + velocitySquared.y * velocitySquared.z)
                    + 20. * velocitySquared.x * velocitySquared.y * velocitySquared.z) / 16.0;

            // Add in the (8) term.
            chi += 81.0 * velocityMagnitudeEight / 128.0
                + 81.0
                    * (velocityEight.x + velocityEight.y + velocityEight.z
                        - (36.0 * velocitySquared.x * velocitySquared.y * velocitySquared.z * velocityMagnitudeSquared
                            + velocityFour.x * velocityFour.y + velocityFour.x * velocityFour.z
                            + velocityFour.y * velocityFour.z)) / 32.0;

            // Multiple whole expression by the density.
            chi *= density;

            util::Vector3D<distribn_t> zeta = util::Vector3D<distribn_t>::Ones() + velocity * 3.0
                + velocitySquared * 9.0 / 2.0 + velocitySquared.PointwiseMultiplication(velocity) * 9.0 / 2.0
                + velocityFour * 27.0 / 8.0;

            zeta.x += CalculateHighOrdersOfZeta<0, 1, 2>(velocity, velocityMagnitudeSquared);
            zeta.y += CalculateHighOrdersOfZeta<1, 2, 0>(velocity, velocityMagnitudeSquared);
            zeta.z += CalculateHighOrdersOfZeta<2, 0, 1>(velocity, velocityMagnitudeSquared);

            for (Direction direction = 0; direction < DmQn::NUMVECTORS; ++direction)
            {
              f_eq[direction] = DmQn::EQMWEIGHTS[direction] * chi
                  * util::NumericalFunctions::IntegerPower(zeta.x, DmQn::CX[direction])
                  * util::NumericalFunctions::IntegerPower(zeta.y, DmQn::CY[direction])
                  * util::NumericalFunctions::IntegerPower(zeta.z, DmQn::CZ[direction]);
            }
          }

          inline static LatticeInfo& GetLatticeInfo()
          {
            if (singletonInfo == NULL)
            {
              util::Vector3D<int> vectors[DmQn::NUMVECTORS];
              Direction inverseVectorIndices[DmQn::NUMVECTORS];

              for (Direction direction = 0; direction < DmQn::NUMVECTORS; ++direction)
              {
                vectors[direction] = util::Vector3D<int>(DmQn::CX[direction], DmQn::CY[direction], DmQn::CZ[direction]);
                inverseVectorIndices[direction] = DmQn::INVERSEDIRECTIONS[direction];
              }

              singletonInfo = new LatticeInfo(DmQn::NUMVECTORS, vectors, inverseVectorIndices);
            }

            return *singletonInfo;
          }

        private:
          inline static distribn_t CalculateStrainRateTensorComponent(const unsigned &iRow,
                                                                      const unsigned &iColumn,
                                                                      const distribn_t &iTau,
                                                                      const distribn_t iFNeq[],
                                                                      const distribn_t &iDensity)
          {
            distribn_t strain_rate_tensor_i_j = 0.0;

            for (Direction vec_index = 0; vec_index < DmQn::NUMVECTORS; vec_index++)
            {
              strain_rate_tensor_i_j +=
                  iFNeq[vec_index]
                      * (DmQn::discreteVelocityVectors[iRow][vec_index]
                          * DmQn::discreteVelocityVectors[iColumn][vec_index]);
            }

            strain_rate_tensor_i_j *= -1.0 / (2.0 * iTau * iDensity * Cs2);

            return strain_rate_tensor_i_j;
          }

          /**
           * Calculate high order of zeta as defined by equation 10 in Chikatamarla et al (PRL, 97, 010201 (2006)
           * @param velocity
           * @param velocityMagnitudeSquared
           * @return
           */
          template<unsigned thisIndex, unsigned otherIndex1, unsigned otherIndex2>
          inline static distribn_t CalculateHighOrdersOfZeta(const util::Vector3D<distribn_t>& velocity,
                                                             distribn_t velocityMagnitudeSquared)
          {
            // Get the velocity components. Note that the naming is to make it easier to follow the
            // paper. ux does not necessarily hold the velocity in the x direction; it's the velocity
            // component in the direction we're calculating zeta for.
            distribn_t ux = velocity[thisIndex], uy = velocity[otherIndex1], uz = velocity[otherIndex2];

            // The 5th order term.
            distribn_t zetaHighOrders = 27.0
                * (util::NumericalFunctions::IntegerPower(ux, 5) - 4. * ux * uy * uy * uz * uz) / 8.0;

            // The 6th order term.
            zetaHighOrders += 81.0 * (util::NumericalFunctions::IntegerPower(ux, 6) - 8. * ux * ux * uy * uy * uz * uz)
                / 16.0;

            // The 7th order term.
            zetaHighOrders += 81.0
                * (util::NumericalFunctions::IntegerPower(ux, 7)
                    + 2. * ux * uy * uy * uz * uz * velocityMagnitudeSquared - 10. * ux * ux * ux * uy * uy * uz * uz)
                / 16.0;

            // The 8th order term.
            zetaHighOrders += 243.0
                * (util::NumericalFunctions::IntegerPower(ux, 8)
                    + 16.0 * ux * ux * uy * uy * uz * uz * (uy * uy + uz * uz)) / 128.0;

            return zetaHighOrders;
          }

          static LatticeInfo* singletonInfo;
      };
    }
  }
}

#endif
