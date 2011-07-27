#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ELBM_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ELBM_H

#include "lb/collisions/implementations/CollisionOperator.h"
#include <cstdlib>

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        class ELBM : public CollisionOperator
        {

          public:
            static void createAlphaArray(const size_t size);

            static void setTau(double* t);

            static void getSiteValues(const distribn_t* f,
                                      distribn_t &density,
                                      distribn_t &v_x,
                                      distribn_t &v_y,
                                      distribn_t &v_z,
                                      distribn_t* f_eq,
                                      const site_t index);

            // WARNING: DOES NOT CALCULATE ALPHA, BECAUSE NON OF THE CURRENT BCS USE IT
            static void getBoundarySiteValues(const distribn_t* f,
                                              const distribn_t &density,
                                              const distribn_t &v_x,
                                              const distribn_t &v_y,
                                              const distribn_t &v_z,
                                              distribn_t* f_eq,
                                              const site_t index);

            static distribn_t getOperatorElement(distribn_t &f_i,
                                                 distribn_t &f_neq_i,
                                                 const LbmParameters* iLbmParams);

          private:
            static double* tau;
            static double* alpha;
            static size_t currentAlphaIndex;
            static double getAlpha(const distribn_t* lFOld,const  distribn_t* lFEq, double prevAlpha);
        };

      }
    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ELBM_H */
