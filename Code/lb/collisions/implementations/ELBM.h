#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ELBM_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ELBM_H

#include "lb/collisions/CollisionOperator.h"

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
            ELBM(const geometry::LatticeData* iLatDat, const lb::LbmParameters* iLbmParams);

            ~ELBM();

            void getSiteValues(const distribn_t* f,
                               distribn_t &density,
                               distribn_t &v_x,
                               distribn_t &v_y,
                               distribn_t &v_z,
                               distribn_t* f_eq,
                               const site_t index);

            // WARNING: DOES NOT CALCULATE ALPHA, BECAUSE NON OF THE CURRENT BCS USE IT
            void getBoundarySiteValues(const distribn_t* f,
                                       const distribn_t &density,
                                       const distribn_t &v_x,
                                       const distribn_t &v_y,
                                       const distribn_t &v_z,
                                       distribn_t* f_eq,
                                       const site_t index);

            distribn_t getOperatorElement(distribn_t &f_i, distribn_t &f_neq_i);

            void Reset(const geometry::LatticeData* iLatDat, const lb::LbmParameters* iLbmParams);

          private:
            double Beta;
            double TwoTau;
            double* alpha;
            size_t currentAlphaIndex;

            double getAlpha(const distribn_t* lFOld, const distribn_t* lFEq, double prevAlpha);
        };

      }
    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ELBM_H */
