#ifndef HEMELB_UNITTESTS_LBTESTS_LBTESTSHELPER_H
#define HEMELB_UNITTESTS_LBTESTS_LBTESTSHELPER_H

#include <cmath>
#include "constants.h"
#include "lb/kernels/BaseKernel.h"
#include "lb/MacroscopicPropertyCache.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      class LbTestsHelper
      {
        public:
          template<typename Lattice>
          static void CalculateRhoVelocity(const distribn_t f[Lattice::NUMVECTORS], distribn_t& rho, distribn_t v[3])
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
          static void CalculateAnsumaliEntropicEqmF(distribn_t density,
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
            lb::HFunction<Lattice> HFunc(f, f_eq);

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
                                    distribn_t expectedFEq[lb::lattices::D3Q15::NUMVECTORS],
                                    std::string id,
                                    lb::kernels::HydroVars<Kernel> &hydroVars,
                                    distribn_t allowedError)
          {
            // Compare density
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Density " + id, expectedDensity, hydroVars.density, allowedError);

            // Compare velocity
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Vx " + id, expectedVx, hydroVars.v_x, allowedError);
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Vy " + id, expectedVy, hydroVars.v_y, allowedError);
            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Vz " + id, expectedVz, hydroVars.v_z, allowedError);

            // Compare equilibrium f
            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              std::stringstream message("FEq ");
              message << ii << " " << id;

              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message.str(),
                                                   expectedFEq[ii],
                                                   hydroVars.GetFEq().f[ii],
                                                   allowedError);
            }
          }

          template<typename LatticeType>
          static void InitialiseAnisotropicTestData(hemelb::geometry::LatticeData* latticeData)
          {
            for (site_t site = 0; site < latticeData->GetLocalFluidSiteCount(); ++site)
            {
              distribn_t* fOld = latticeData->GetSite(site).GetFOld<LatticeType>();
              InitialiseAnisotropicTestData<LatticeType>(site, fOld);
            }
          }

          template<typename LatticeType>
          static void InitialiseAnisotropicTestData(site_t site, distribn_t* distribution)
          {
            for (unsigned int direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
            {
              distribution[direction] = ((distribn_t) (direction + 1)) / 10.0 + ((distribn_t) (site)) / 100.0;
            }
          }

          template<class Kernel>
          static void CalculateRegularisedCollision(geometry::LatticeData* latDat,
                                                    lb::LbmParameters* lbmParams,
                                                    site_t siteIndex,
                                                    lb::kernels::HydroVars<Kernel>& hydroVars,
                                                    distribn_t* const fPostCollision)
          {
            const distribn_t * const fPreCollision = hydroVars.f;

            // To evaluate PI, first let unknown particle populations take value given by bounce-back of off-equilibrium parts
            // (fi = fiEq + fopp(i) - fopp(i)Eq)
            distribn_t fTemp[lb::lattices::D3Q15::NUMVECTORS];

            for (unsigned l = 0; l < lb::lattices::D3Q15::NUMVECTORS; ++l)
            {
              if (latDat->GetSite(siteIndex).HasBoundary(lb::lattices::D3Q15::INVERSEDIRECTIONS[l]))
              {
                fTemp[l] = hydroVars.GetFEq().f[l] + fPreCollision[lb::lattices::D3Q15::INVERSEDIRECTIONS[l]]
                    - hydroVars.GetFEq().f[lb::lattices::D3Q15::INVERSEDIRECTIONS[l]];
              }
              else
              {
                fTemp[l] = fPreCollision[l];
              }
            }

            distribn_t f_neq[lb::lattices::D3Q15::NUMVECTORS];
            for (unsigned l = 0; l < lb::lattices::D3Q15::NUMVECTORS; ++l)
            {
              f_neq[l] = fTemp[l] - hydroVars.GetFEq().f[l];
            }

            // Pi = sum_i e_i e_i f_i
            // zeta = Pi / 2 (Cs^4)
            Order2Tensor zeta = lb::lattices::D3Q15::CalculatePiTensor(f_neq);

            for (int m = 0; m < 3; m++)
            {
              for (int n = 0; n < 3; n++)
              {
                zeta[m][n] /= (2.0 * Cs2 * Cs2);
              }
            }

            // chi = Cs^2 I : zeta
            const distribn_t chi = Cs2 * (zeta[0][0] + zeta[1][1] + zeta[2][2]);

            const int *Cs[3] = { lb::lattices::D3Q15::CX, lb::lattices::D3Q15::CY, lb::lattices::D3Q15::CZ };

            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              // According to Latt & Chopard (Physical Review E77, 2008),
              // f_neq[i] = (LatticeWeight[i] / (2 Cs^4)) *
              //            Q_i : Pi(n_eq)
              // Where Q_i = c_i c_i - Cs^2 I
              // and Pi(n_eq) = Sum{i} (c_i c_i f_i)
              //
              // We pre-compute zeta = Pi(neq) / (2 Cs^4)
              //             and chi =  Cs^2 I : zeta
              // Hence we can compute f_neq[i] = LatticeWeight[i] * ((c_i c_i) : zeta - chi)
              f_neq[ii] = -chi;

              for (int aa = 0; aa < 3; ++aa)
              {
                for (int bb = 0; bb < 3; ++bb)
                {
                  f_neq[ii] += (float(Cs[aa][ii] * Cs[bb][ii])) * zeta[aa][bb];
                }
              }

              f_neq[ii] *= lb::lattices::D3Q15::EQMWEIGHTS[ii];

              /*
               * Newly constructed distribution function:
               *    g_i = f^{eq}_i + f^{neq}_i
               *
               * Collision step:
               *    f^{+}_i = g_i + w (g_i - f^{eq}_i)
               *            = f^{eq}_i + (1+w) f^{neq}_i
               */
              fPostCollision[ii] = hydroVars.GetFEq().f[ii] + (1.0 + lbmParams->GetOmega()) * f_neq[ii];
            }

          }

          /**
           * Updates a property cache for the macroscopic variables selected. This should have
           * identical behaviour to in UpdateSiteMinsAndMaxes in BaseStreamer where the cache
           * is normally populated.
           * @param latDat
           * @param cache
           * @param simState
           */
          template<typename Lattice>
          static void UpdatePropertyCache(geometry::LatticeData& latDat,
                                          lb::MacroscopicPropertyCache& cache,
                                          lb::SimulationState& simState)
          {
            for (site_t site = 0; site < latDat.GetLocalFluidSiteCount(); ++site)
            {
              distribn_t density, feq[Lattice::NUMVECTORS];
              util::Vector3D<distribn_t> velocity;

              Lattice::CalculateDensityVelocityFEq(latDat.GetSite(site).GetFOld<Lattice>(),
                                                   density,
                                                   velocity[0],
                                                   velocity[1],
                                                   velocity[2],
                                                   feq);

              if (cache.densityCache.RequiresRefresh())
              {
                cache.densityCache.Put(site, density);
              }
              if (cache.velocityCache.RequiresRefresh())
              {
                cache.velocityCache.Put(site, velocity / density);
              }

              // TODO stress cache filling not yet implemented.
            }
          }
      };
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_LBTESTSHELPER_H */
