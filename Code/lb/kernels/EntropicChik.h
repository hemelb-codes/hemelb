#ifndef HEMELB_LB_KERNELS_ENTROPICCHIK_H
#define HEMELB_LB_KERNELS_ENTROPICCHIK_H

#include <cstdlib>

#include "lb/kernels/Entropic.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      // We have to declare this up here in order for it to be used as a template parameter in the
      // following declaration. Moving the template specialisation to the bottom of the file would
      // prevent it from being used as the HydroVars for this kernel.
      template<class LatticeType> class EntropicChik;

      template<class LatticeType>
      struct HydroVars<EntropicChik<LatticeType> > : public HydroVarsBase<LatticeType>
      {
        public:
          HydroVars(const distribn_t* const f) :
              HydroVarsBase<LatticeType>(f)
          {

          }

          site_t index;
      };

      /**
       * EntropicChik: This class implements the entropic kernel, as per Chitakamarla et al.
       */
      template<class LatticeType>
      class EntropicChik : public BaseKernel<EntropicChik<LatticeType>, LatticeType>
                           , public Entropic<LatticeType>
      {
        public:
          /**
           * Constructor, passes parameters onto the base class.
           * @param initParams
           */
          EntropicChik(InitParams& initParams) :
              Entropic<LatticeType>(&initParams)
          {
          }

          /**
           * Calculates the density and velocity for the given f. Then calculates the
           * equilibrium distribution as described by Chikatamarla.
           * @param hydroVars
           * @param index
           */
          inline void DoCalculateDensityVelocityFeq(HydroVars<EntropicChik<LatticeType> >& hydroVars, site_t index)
          {
            hydroVars.index = index;
            LatticeType::CalculateDensityAndVelocity(hydroVars.f,
                                                     hydroVars.density,
                                                     hydroVars.v_x,
                                                     hydroVars.v_y,
                                                     hydroVars.v_z);
            LatticeType::CalculateEntropicFeqChik(hydroVars.density,
                                                  hydroVars.v_x,
                                                  hydroVars.v_y,
                                                  hydroVars.v_z,
                                                  hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
          }

          /**
           * Calculates the equilibrium f distribution for the given density and velocity, as
           * described by Chikatamarla.
           * @param hydroVars
           * @param index
           */
          inline void DoCalculateFeq(HydroVars<EntropicChik<LatticeType> >& hydroVars, site_t index)
          {
            hydroVars.index = index;
            LatticeType::CalculateEntropicFeqChik(hydroVars.density,
                                                  hydroVars.v_x,
                                                  hydroVars.v_y,
                                                  hydroVars.v_z,
                                                  hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
          }
      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_ENTROPICANSUMALI_H */
