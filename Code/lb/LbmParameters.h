#ifndef HEMELB_LB_LBMPARAMETERS_H
#define HEMELB_LB_LBMPARAMETERS_H

#include <cmath>
#include "constants.h"
#include <vector>
#include <cassert>
#include "D3Q15.h"
#include "lb/kernels/momentBasis/DHumieresD3Q15MRTBasis.h"

namespace hemelb
{
  namespace lb
  {
    enum StressTypes
    {
      VonMises = 0,
      ShearStress = 1,
      IgnoreStress = 2
    };

    struct LbmParameters
    {
      public:
        LbmParameters(distribn_t timeStepLength, distribn_t voxelSize)
        {
          Update(timeStepLength, voxelSize);
        }

        void Update(distribn_t timeStepLength, distribn_t voxelSize)
        {
          timestep = timeStepLength;
          tau = 0.5
              + (timeStepLength * BLOOD_VISCOSITY_Pa_s / BLOOD_DENSITY_Kg_per_m3)
                  / (Cs2 * voxelSize * voxelSize);

          omega = -1.0 / tau;
          stressParameter = (1.0 - 1.0 / (2.0 * tau)) / sqrt(2.0);
          beta = -1.0 / (2.0 * tau);

          // Relaxation values taken from d'Humieres 2002, except for the kinematic viscosity where the usual tau formula is used.
          mrtRelaxationParameters.clear();
          mrtRelaxationParameters.push_back(1.6); // e (s1)
          mrtRelaxationParameters.push_back(1.2); // epsilon (s2)
          mrtRelaxationParameters.push_back(1.6); // q_x (s4)
          mrtRelaxationParameters.push_back(1.6); // q_y (s4)
          mrtRelaxationParameters.push_back(1.6); // q_z (s4)
          mrtRelaxationParameters.push_back(1.0 / tau); // 3p_xx (s9)
          mrtRelaxationParameters.push_back(1.0 / tau); // p_ww (s9)
          mrtRelaxationParameters.push_back(1.0 / tau); // p_xy (s11 = s9)
          mrtRelaxationParameters.push_back(1.0 / tau); // p_yz (s11 = s9)
          mrtRelaxationParameters.push_back(1.0 / tau); // p_zx (s11 = s9)
          mrtRelaxationParameters.push_back(1.2); // m_xyz (s14)
          assert(mrtRelaxationParameters.size() == kernels::momentBasis::DHumieresD3Q15MRTBasis::NUM_KINETIC_MOMENTS);
        }

        distribn_t GetTimeStep() const
        {
          return timestep;
        }

        distribn_t GetOmega() const
        {
          return omega;
        }

        distribn_t GetTau() const
        {
          return tau;
        }

        distribn_t GetStressParameter() const
        {
          return stressParameter;
        }

        distribn_t GetBeta() const
        {
          return beta;
        }

        const std::vector<distribn_t>& GetMrtRelaxationParameters() const
        {
          return mrtRelaxationParameters;
        }

        void SetMrtRelaxationParameters(std::vector<distribn_t>& newRelaxationParameters)
        {
          mrtRelaxationParameters = newRelaxationParameters;
        }

        StressTypes StressType;

      private:
        distribn_t timestep;
        distribn_t omega;
        distribn_t tau;
        distribn_t stressParameter;
        distribn_t beta; ///< Viscous dissipation in ELBM
        std::vector<distribn_t> mrtRelaxationParameters; ///< Relaxation paremeters used by the MRT kernel.
    };
  }
}

#endif //HEMELB_LB_LBMPARAMETERS_H
