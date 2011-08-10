#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_LBGKNN_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_LBGKNN_H

#include "lb/collisions/implementations/LBGK.h"

#include "lb/rheology_models/CassonRheologyModel.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        /*
         * Class extending the original BGK collision operator to support non-newtonian
         * fluids. Overrides the relevant methods and implements support for relaxation
         * time not constant across the domain.
         */
        template<class tNonNewtonianModel>
        class LBGKNN : public LBGK
        {
          public:

            /*
             * Creates the relevant data structures for allowing different values of tau
             * across the domain.
             */
            static void createTauArray(const site_t& iSize, const double& iDefaultTau);

            static void getSiteValues(const distribn_t* f,
                                      distribn_t &density,
                                      distribn_t &v_x,
                                      distribn_t &v_y,
                                      distribn_t &v_z,
                                      distribn_t* f_eq,
                                      const site_t index);

//            static void getBoundarySiteValues(const distribn_t* f,
//                                              const distribn_t &density,
//                                              const distribn_t &v_x,
//                                              const distribn_t &v_y,
//                                              const distribn_t &v_z,
//                                              distribn_t* f_eq,
//                                              const site_t index);
//
//            static void doPostCalculations(const distribn_t* f,
//                                           const geometry::LatticeData* bLatDat,
//                                           const site_t inde);

            static distribn_t getOperatorElement(distribn_t &f_i,
                                                 distribn_t &f_neq_i,
                                                 const LbmParameters* iLbmParams);

          private:
            /* Vector containing the current value of tau for each site in the domain */
            static std::vector<double> mTau;

            /*
             * Since getOperatorElement doesn't get the index of the site it is working with,
             * or the local density, we store them when getSiteValues is called (which is
             * supposed to happen before)
             */
            static size_t mCurrentTauIndex;
            static distribn_t mCurrentDensity;

        };
      }
    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_LBGKNN_H */
