#ifndef HEMELB_UNITTESTS_LBTESTS_KERNELTESTSHELPER_H
#define HEMELB_UNITTESTS_LBTESTS_KERNELTESTSHELPER_H

#include <cmath>
#include "constants.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      class KernelTestsHelper
      {
        public:
          template<typename Lattice>
          static void CalculateRhoVelocity(const distribn_t f[Lattice::NUMVECTORS],
                                           distribn_t& rho,
                                           distribn_t v[3])
          {
            rho = 0.0;

            v[0] = v[1] = v[2] = 0.0;
            for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
            {
              rho += f[ii];
              v[0] += f[ii] * (float) Lattice::CX[ii];
              v[1] += f[ii] * (float) Lattice::CY[ii];
              v[2] += f[ii] * (float) Lattice::CZ[ii];
            }
          }

          template<typename Lattice>
          static void CalculateVelocity(const distribn_t f[Lattice::NUMVECTORS], distribn_t v[3])
          {
            v[0] = v[1] = v[2] = 0.0;
            for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
            {
              v[0] += f[ii] * (float) Lattice::CX[ii];
              v[1] += f[ii] * (float) Lattice::CY[ii];
              v[2] += f[ii] * (float) Lattice::CZ[ii];
            }
          }

          template<typename Lattice>
          static void CalculateEntropicEqmF(distribn_t density,
                                            distribn_t v_x,
                                            distribn_t v_y,
                                            distribn_t v_z,
                                            distribn_t f[Lattice::NUMVECTORS])
          {
            // Calculate velocity.
            distribn_t u[3] = { v_x / density, v_y / density, v_z / density };

            distribn_t B[3];
            distribn_t term1[3];
            distribn_t term2[3];
            for (int ii = 0; ii < 3; ++ii)
            {
              B[ii] = sqrt(1.0 + 3.0 * u[ii] * u[ii]);
              term1[ii] = 2.0 - B[ii];
              term2[ii] = (2.0 * u[ii] + B[ii]) / (1.0 - u[ii]);
            }

            for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
            {
              f[ii] = density * Lattice::EQMWEIGHTS[ii];

              f[ii] *= term1[0] * pow(term2[0], (double) Lattice::CX[ii]);
              f[ii] *= term1[1] * pow(term2[1], (double) Lattice::CY[ii]);
              f[ii] *= term1[2] * pow(term2[2], (double) Lattice::CZ[ii]);
            }
          }

          template<typename Lattice>
          static void CalculateLBGKEqmF(distribn_t density,
                                        distribn_t v_x,
                                        distribn_t v_y,
                                        distribn_t v_z,
                                        distribn_t f[Lattice::NUMVECTORS])
          {
            for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
            {
              // Calculate the dot-product of the velocity with the direction vector.
              distribn_t vSum = v_x * (float) Lattice::CX[ii] + v_y * (float) Lattice::CY[ii]
                  + v_z * (float) Lattice::CZ[ii];

              // Calculate the squared magnitude of the velocity.
              distribn_t v2Sum = v_x * v_x + v_y * v_y + v_z * v_z;

              // F eqm = density proportional component...
              f[ii] = density;

              // ... - v^2 component...
              f[ii] -= ( (3.0 / 2.0) * v2Sum / density);

              // ... + v^1 component
              f[ii] += 3.0 * vSum + (9.0 / 2.0) * vSum * vSum / density;

              // Multiply by eqm weight.
              f[ii] *= Lattice::EQMWEIGHTS[ii];
            }
          }

          template<typename Lattice>
          static void CalculateEntropicCollision(const distribn_t f[Lattice::NUMVECTORS],
                                                 const distribn_t f_eq[Lattice::NUMVECTORS],
                                                 distribn_t tau,
                                                 distribn_t beta,
                                                 distribn_t f_collided[Lattice::NUMVECTORS])
          {
            lb::HFunction HFunc(f, f_eq);

            distribn_t alpha = hemelb::util::NumericalMethods::NewtonRaphson(&HFunc, 2.0, 1.0E-100);

            for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
            {
              f_collided[ii] = f[ii] + alpha * beta * (f[ii] - f_eq[ii]);
            }
          }

          template<typename Lattice>
          static void CalculateLBGKCollision(const distribn_t f[Lattice::NUMVECTORS],
                                             const distribn_t f_eq[Lattice::NUMVECTORS],
                                             distribn_t omega,
                                             distribn_t f_collided[Lattice::NUMVECTORS])
          {
            for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
            {
              f_collided[ii] = f[ii] + omega * (f[ii] - f_eq[ii]);
            }
          }

          template<typename Kernel>
          static void CompareHydros(distribn_t expectedDensity,
                                    distribn_t expectedVx,
                                    distribn_t expectedVy,
                                    distribn_t expectedVz,
                                    distribn_t expectedFEq[D3Q15::NUMVECTORS],
                                    std::string id,
                                    lb::kernels::HydroVars<Kernel> &hydroVars,
                                    distribn_t allowedError)
          {
            // Compare density
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Density " + id,
                                                 expectedDensity,
                                                 hydroVars.density,
                                                 allowedError);

            // Compare velocity
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Vx " + id,
                                                 expectedVx,
                                                 hydroVars.v_x,
                                                 allowedError);
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Vy " + id,
                                                 expectedVy,
                                                 hydroVars.v_y,
                                                 allowedError);
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Vz " + id,
                                                 expectedVz,
                                                 hydroVars.v_z,
                                                 allowedError);

            // Compare equilibrium f
            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              std::stringstream message("FEq ");
              message << ii << " " << id;

              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message.str(),
                                                   expectedFEq[ii],
                                                   hydroVars.f_eq[ii],
                                                   allowedError);
            }
          }
      };
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_KERNELTESTSHELPER_H */
