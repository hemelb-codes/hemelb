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
#ifdef HEMELB_USE_SSE3
  #include <immintrin.h>
#endif
#include <mpi.h>
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
          
         #ifdef HEMELB_USE_SSE3
          /**
           * Calculates density and momentum using SSE3 intrinsics.
           * If the lattice has an odd number of vectors (directions), 
           * the last element is processed first 
           * 
           * The reductions are calculated in two streams (the loop is virtually 
           * twice unrolled), followed by horizontal adds to sum two partial results together
           * 
           * @param f
           * @param density
           * @param momentum_x
           * @param momentum_y
           * @param momentum_z
           */
          inline static void CalculateDensityAndMomentum(const distribn_t f[],
                                                       distribn_t &density,
                                                       distribn_t &momentum_x,
                                                       distribn_t &momentum_y,
                                                       distribn_t &momentum_z)
          {                        
            // SSE2 accumulator registers containing a pair of double values
            __m128d density_SSE2;
            __m128d momentum_x_SSE2;
            __m128d momentum_y_SSE2;
            __m128d momentum_z_SSE2;

            // set the loop boundary to the highest even number  =< DmQn::NUMVECTORS
            Direction numVect2 = ((DmQn::NUMVECTORS >> 1) << 1); 

            // process the 15/19/27th element first
            if (DmQn::NUMVECTORS != numVect2) 
            {
              // the first double is set to the result of the last element, the second double to zero
              density_SSE2 = _mm_set_pd(f[DmQn::NUMVECTORS - 1], 0.0);
              momentum_x_SSE2 = _mm_set_pd(DmQn::CXD[DmQn::NUMVECTORS - 1] * f[DmQn::NUMVECTORS - 1], 0.0);
              momentum_y_SSE2 = _mm_set_pd(DmQn::CYD[DmQn::NUMVECTORS - 1] * f[DmQn::NUMVECTORS - 1], 0.0);
              momentum_z_SSE2 = _mm_set_pd(DmQn::CZD[DmQn::NUMVECTORS - 1] * f[DmQn::NUMVECTORS - 1], 0.0);
            }
            else
            {
              // set the SSE accumulators to zero
              density_SSE2 = _mm_set1_pd(0.0);
              momentum_x_SSE2 = _mm_set1_pd(0.0);
              momentum_y_SSE2 = _mm_set1_pd(0.0);
              momentum_z_SSE2 = _mm_set1_pd(0.0);
            }

            //SSE loop processing two elements at once
            for (Direction direction = 0; direction < numVect2; direction += 2)
            {

              // f is not aligned, loadu has to be used,
              // CXD, CYD, CZD are supposed to be 16B aligned
              const __m128d f_SSE2 = _mm_loadu_pd(&f[direction]);
              const __m128d CX_SSE2 = _mm_load_pd(&DmQn::CXD[direction]);
              const __m128d CY_SSE2 = _mm_load_pd(&DmQn::CYD[direction]);
              const __m128d CZ_SSE2 = _mm_load_pd(&DmQn::CZD[direction]);

              // density += f[i]
              density_SSE2 = _mm_add_pd(density_SSE2, f_SSE2);

              // momentum_x += CX[i] * f[i]]
              momentum_x_SSE2 = _mm_add_pd(momentum_x_SSE2, _mm_mul_pd(CX_SSE2, f_SSE2));
              momentum_y_SSE2 = _mm_add_pd(momentum_y_SSE2, _mm_mul_pd(CY_SSE2, f_SSE2));
              momentum_z_SSE2 = _mm_add_pd(momentum_z_SSE2, _mm_mul_pd(CZ_SSE2, f_SSE2));
            }

            // horizontal adds to sum partial results and store them back to memory
            _mm_store_sd(&density,    _mm_hadd_pd(density_SSE2, density_SSE2));
            _mm_store_sd(&momentum_x, _mm_hadd_pd(momentum_x_SSE2, momentum_x_SSE2));
            _mm_store_sd(&momentum_y, _mm_hadd_pd(momentum_y_SSE2, momentum_y_SSE2));
            _mm_store_sd(&momentum_z, _mm_hadd_pd(momentum_z_SSE2, momentum_z_SSE2));

          }
          
         #else

          /**
           * Calculates density and momentum, the original non-SSE version
           * @param f
           * @param density
           * @param momentum_x
           * @param momentum_y
           * @param momentum_z
           */
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
          #endif                   

          #ifdef HEMELB_USE_SSE3
          /**           
           * Calculates Feq using SSE3 intrinsics.
           * If the lattice has an odd number of vectors (directions), 
           * the last element is processed using scalar arithmetics
           * 
           * The reductions are calculated in two streams (the loop is virtually 
           * twice unrolled). Some invariants are merged together
           * 
           * @param density
           * @param momentum_x
           * @param momentum_y
           * @param momentum_z
           * @param f_eq
           */
          inline static void CalculateFeq(const distribn_t &density,
                                          const distribn_t &momentum_x,
                                          const distribn_t &momentum_y,
                                          const distribn_t &momentum_z,
                                          distribn_t f_eq[])
          {
                       
            // merge some constants and invariants and populate SSE registers by them                        
            const distribn_t threeHalvesOfMomentumMagnitudeSquared = (3./2.) * (momentum_x * momentum_x + momentum_y * momentum_y
                  + momentum_z * momentum_z);            
            const __m128d threeHalvesOfMomentumMagnitudeSquared_SSE2 = _mm_set1_pd(threeHalvesOfMomentumMagnitudeSquared);            
            
            
            const distribn_t density_1 = 1. / density;            
            const __m128d density_1_SSE2 = _mm_set1_pd(density_1);
            const __m128d density_SSE2   = _mm_set1_pd(density);
            
            const __m128d momentum_x_SSE2 = _mm_set1_pd(momentum_x);
            const __m128d momentum_y_SSE2 = _mm_set1_pd(momentum_y);
            const __m128d momentum_z_SSE2 = _mm_set1_pd(momentum_z);
            
            const distribn_t nineHalvesOfDensity_1 =  (9. / 2.) * density_1;
            const __m128d nineOnTwoDensity_1_SSE2  = _mm_set1_pd(nineHalvesOfDensity_1) ;
            const __m128d three_SSE2      = _mm_set1_pd(3.) ;
            
            
            // sse loop (the loop is virtually twice unrolled)
            Direction numVect2 = ((DmQn::NUMVECTORS >> 1) << 1); 
            for (Direction i = 0; i < numVect2 ; i+=2)            
            {
              // mom_dot_ei = DmQn::CX[i] * momentum_x + DmQn::CY[i] * momentum_y + DmQn::CZ[i] * momentum_z;
              const __m128d CXD_momentum_x_SSE2 = _mm_mul_pd(_mm_load_pd(&DmQn::CXD[i]),momentum_x_SSE2);            
              const __m128d CYD_momentum_y_SSE2 = _mm_mul_pd(_mm_load_pd(&DmQn::CYD[i]),momentum_y_SSE2);            
              const __m128d CZD_momentum_z_SSE2 = _mm_mul_pd(_mm_load_pd(&DmQn::CZD[i]),momentum_z_SSE2);            
              
              const __m128d EQMWEIGHTS_SSE2 = _mm_load_pd(&DmQn::EQMWEIGHTS[i]);
              
              const __m128d mom_dot_ei_SSE2 =  _mm_add_pd(
                                                    _mm_add_pd(CXD_momentum_x_SSE2, CYD_momentum_y_SSE2), 
                                                    CZD_momentum_z_SSE2
                                               );
                                                
              
              //  (density - (3. / 2.) * momentumMagnitudeSquared * density_1
              const __m128d tmp1 = _mm_sub_pd(density_SSE2,
                                              _mm_mul_pd(threeHalvesOfMomentumMagnitudeSquared_SSE2, density_1_SSE2 )
                                             );
              
              // (9. / 2.) * density_1 * mom_dot_ei * mom_dot_ei
              const __m128d tmp2 = (_mm_mul_pd(
                                              nineOnTwoDensity_1_SSE2,
                                              _mm_mul_pd(mom_dot_ei_SSE2, mom_dot_ei_SSE2)
                                              )
                                   );
              // 3. * mom_dot_ei);
              const __m128d tmp3 = _mm_mul_pd(three_SSE2, mom_dot_ei_SSE2);
              
              
              
              __m128d tmp4 = _mm_add_pd(tmp1, tmp2);
                      tmp4 = _mm_add_pd(tmp4, tmp3);
                                          
              // f_eq is not 16B aligned
              _mm_storeu_pd(&f_eq[i],_mm_mul_pd(EQMWEIGHTS_SSE2,tmp4));
            }
              
            
            // do the odd element (15/19/27) 
            if (DmQn::NUMVECTORS != numVect2) // constants are reduced
            {
                          
              const distribn_t mom_dot_ei = DmQn::CX[DmQn::NUMVECTORS-1] * momentum_x
                  + DmQn::CY[DmQn::NUMVECTORS-1] * momentum_y + DmQn::CZ[DmQn::NUMVECTORS-1] * momentum_z;

              f_eq[DmQn::NUMVECTORS-1] = DmQn::EQMWEIGHTS[DmQn::NUMVECTORS - 1]
                  * (density -  threeHalvesOfMomentumMagnitudeSquared  * density_1
                      + nineHalvesOfDensity_1 * (mom_dot_ei * mom_dot_ei)
                      + 3. * mom_dot_ei);
             
            }
          }
         #else                    
          
          /**
           * Calculate Feq, the orginal version
           * @param density
           * @param momentum_x
           * @param momentum_y
           * @param momentum_z
           * @param f_eq
           */
          inline static void CalculateFeq(const distribn_t &density,
                                          const distribn_t &momentum_x,
                                          const distribn_t &momentum_y,
                                          const distribn_t &momentum_z,
                                          distribn_t f_eq[])
          {
            const distribn_t density_1 = 1. / density;
            const distribn_t momentumMagnitudeSquared = momentum_x * momentum_x + momentum_y * momentum_y
                + momentum_z * momentum_z;

            for (Direction i = 0; i < DmQn::NUMVECTORS; ++i)
            {
              const distribn_t mom_dot_ei = DmQn::CX[i] * momentum_x
                  + DmQn::CY[i] * momentum_y + DmQn::CZ[i] * momentum_z;

              f_eq[i] = DmQn::EQMWEIGHTS[i]
                  * (density - (3. / 2.) * momentumMagnitudeSquared * density_1
                      + (9. / 2.) * density_1 * mom_dot_ei * mom_dot_ei
                      + 3. * mom_dot_ei);
            }
          }
          #endif
          
                                        
          // Calculate density, momentum and the equilibrium distribution
          // functions according to the D3Q15 model.  The calculated momentum_x, momentum_y
          // and momentum_z are actually density * velocity, because we are using the
          // compressible model.
          inline static void CalculateDensityMomentumFEq(const distribn_t f[],
                                                         distribn_t &density,
                                                         distribn_t &momentum_x,
                                                         distribn_t &momentum_y,
                                                         distribn_t &momentum_z,
                                                         distribn_t &velocity_x,
                                                         distribn_t &velocity_y,
                                                         distribn_t &velocity_z,
                                                         distribn_t f_eq[])
          {                       
            CalculateDensityAndMomentum(f, density, momentum_x, momentum_y, momentum_z);            
                       
            velocity_x = momentum_x / density;
            velocity_y = momentum_y / density;
            velocity_z = momentum_z / density;
            
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
           * @param traction traction vector at a given point
           */
          inline static void CalculateTractionOnAPoint(const distribn_t density,
                                                       const distribn_t tau,
                                                       const distribn_t fNonEquilibrium[],
                                                       const util::Vector3D<Dimensionless>& wallNormal,
                                                       util::Vector3D<LatticeStress>& traction)
          {
            util::Matrix3D sigma;
            CalculateStressTensor(density, tau, fNonEquilibrium, sigma);

            // Multiply the stress tensor by the surface normal
            sigma.timesVector(wallNormal, traction);
          }

          /**
           * Calculates the projection of the traction vector on the plane tangential to the geometry surface (defined by current surface
           * point and the normal vector provided). This is done with the following formula:
           *
           *    \vec{t_tan} = \vec{t} - dot(\vec{t}, \vec{n})*\vec{n}
           *
           * where t is the traction vector (see CalculateTractionOnAPoint for definition) and n is the normal
           *
           * @param density density at a given site
           * @param tau relaxation time
           * @param fNonEquilibrium non equilibrium part of the distribution function
           * @param wallNormal wall normal at a given point
           * @param tractionTangentialComponent tangential projection of the traction vector
           */
          inline static void CalculateTangentialProjectionTraction(const distribn_t density,
                                                                   const distribn_t tau,
                                                                   const distribn_t fNonEquilibrium[],
                                                                   const util::Vector3D<Dimensionless>& wallNormal,
                                                                   util::Vector3D<LatticeStress>& tractionTangentialComponent)
          {
            util::Vector3D<LatticeStress> traction;
            CalculateTractionOnAPoint(density, tau, fNonEquilibrium, wallNormal, traction);

            LatticeStress magnitudeNormalProjectionTraction = traction.Dot(wallNormal);

            tractionTangentialComponent = traction - wallNormal * magnitudeNormalProjectionTraction;
          }

          /**
           * Calculate the full stress tensor at a given fluid site (including both pressure and deviatoric part)
           *
           * The stress tensor is assembled based on the formula:
           *
           *    \sigma = p*I + 2*\mu*S = p*I - \Pi^{(neq)}
           *
           * where p is hydrodynamic pressure, I is the identity tensor, S is the strain rate tensor, and \mu is the
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

            // Add the pressure component to the stress tensor. The reference pressure given
            // by the REFERENCE_PRESSURE_mmHg constant is mapped to rho=1. Here we subtract 1
            // and when the tensor is turned into physical units REFERENCE_PRESSURE_mmHg will
            // be added.
            LatticePressure pressure = (density - 1) * Cs2;
            stressTensor.addDiagonal(pressure);
          }

          /**
           * The magnitude of the tangential component of the shear stress acting on the
           * wall (i.e. tangential component of the traction vector). For this method to
           * make sense f has to be the non equilibrium part of a distribution function.
           *
           * The stress tensor computed in this method only includes the deviatoric part
           * and not the component corresponding to the pressure. Do not use the intermediate
           * traction value unless you understand the implications (use CalculateTractionOnAPoint
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

          inline static bool IsLatticeCompressible()
          {
            return true;
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
