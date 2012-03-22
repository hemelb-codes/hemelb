#ifndef HEMELB_LB_KERNELS_MRT_H
#define HEMELB_LB_KERNELS_MRT_H

#include "lb/kernels/BaseKernel.h"
#include "lb/SimulationState.h"
#include <cassert>
#include <cmath>

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      // Forward declaration needed by the struct
      template<class MomentBasis> class MRT;

      template<class MomentBasis>
      struct HydroVars<MRT<MomentBasis> > : public HydroVarsBase<typename MomentBasis::Lattice>
      {
        public:
          HydroVars(const distribn_t* const f) :
              HydroVarsBase<typename MomentBasis::Lattice>(f)
          {
          }

          /** Equilibrium velocity distribution in the momentum space. */
          distribn_t m_neq[MomentBasis::NUM_KINETIC_MOMENTS];
      };

      /**
       * This class implements the Multiple Relaxation Time (MRT) collision operator.
       *
       *  \Omega(f) = - M^{-1} * \hat{S} * M (f - f_{eq})
       *            = - M^T * (M * M^T)^{-1} * \hat{S} * m_{neq}
       *
       *  where f is the velocity distribution function, \Omega() is the collision operator,
       *  M is the momentum space basis, m=Mf is the vector of momentums, \hat{S} is the
       *  collision matrix, and {m,f}_{eq} and {m,f}_{neq} are the equilibrium and non-equilibrium
       *  versions of {m,f}.
       *
       *  (M * M^T)^{-1} and \hat{S} are diagonal matrices.
       */
      template<class MomentBasis>
      class MRT : public BaseKernel<MRT<MomentBasis>, typename MomentBasis::Lattice>
      {
        public:

          MRT(InitParams& initParams)
          {
            InitState(initParams);
          }

          inline void DoCalculateDensityVelocityFeq(HydroVars<MRT>& hydroVars, site_t index)
          {
            MomentBasis::Lattice::CalculateDensityVelocityFEq(hydroVars.f,
                                                              hydroVars.density,
                                                              hydroVars.v_x,
                                                              hydroVars.v_y,
                                                              hydroVars.v_z,
                                                              hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < MomentBasis::Lattice::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }

            /** @todo #61 consider computing m_neq directly in the momentum space. See d'Humieres 2002. */
            MomentBasis::ProjectVelsIntoMomentSpace(hydroVars.f_neq.f, hydroVars.m_neq);
          }

          inline void DoCalculateFeq(HydroVars<MRT>& hydroVars, site_t index)
          {
            MomentBasis::Lattice::CalculateFeq(hydroVars.density,
                                               hydroVars.v_x,
                                               hydroVars.v_y,
                                               hydroVars.v_z,
                                               hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < MomentBasis::Lattice::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }

            /** @todo #61 consider computing m_neq directly in the momentum space. See d'Humieres 2002. */
            MomentBasis::ProjectVelsIntoMomentSpace(hydroVars.f_neq.f, hydroVars.m_neq);
          }

          inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<MRT>& hydroVars)
          {
            for (Direction direction = 0; direction < MomentBasis::Lattice::NUMVECTORS; ++direction)
            {
              /** @todo #61 many optimisations possible (and necessary!).
               *  - Store the product of REDUCED_MOMENT_BASIS and 1/BASIS_TIMES_BASIS_TRANSPOSED instead of REDUCED_MOMENT_BASIS
               *  - Compute the loop below as a matrix product in DoCalculate*, alternatively we could consider reimplementing DoCollide to work with whole arrays (consider libraries boost::ublas or Armadillo)
               */
              distribn_t collision = 0.;
              for (unsigned momentIndex = 0; momentIndex < MomentBasis::NUM_KINETIC_MOMENTS; momentIndex++)
              {
                collision += (collisionMatrix[momentIndex] / MomentBasis::BASIS_TIMES_BASIS_TRANSPOSED[momentIndex])
                    * MomentBasis::REDUCED_MOMENT_BASIS[momentIndex][direction] * hydroVars.m_neq[momentIndex];
              }
              hydroVars.SetFPostCollision(direction, hydroVars.f[direction] - collision);
            }
          }

          inline void DoReset(InitParams* initParams)
          {
            InitState(*initParams);
          }

          /**
           *  This method is used in unit testing in order to make an MRT kernel behave as LBGK, regardless of the
           *  moment basis, by setting all the relaxation parameters to be the same.
           */
          void SetMrtRelaxationParameters(std::vector<distribn_t>& newRelaxationParameters)
          {
            assert(newRelaxationParameters.size() == MomentBasis::NUM_KINETIC_MOMENTS);
            collisionMatrix = newRelaxationParameters;
          }

        private:
          /** MRT collision matrix (\hat{S}, diagonal). It corresponds to the inverse of the relaxation time for each mode. */
          std::vector<distribn_t> collisionMatrix;

          /**
           *  Helper method to set/update member variables. Called from the constructor and Reset()
           *
           *  @param initParams struct used to store variables required for initialisation of various operators
           */
          void InitState(const InitParams& initParams)
          {
            MomentBasis::SetUpCollisionMatrix(collisionMatrix, initParams.lbmParams->GetTau());
          }

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_MRT_H */
