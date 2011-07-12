#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ELBM_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ELBM_H

#include "topology/NetworkTopology.h"
#include "lb/collisions/implementations/CollisionOperator.h"
#include "lb/collisions/implementations/HFunction.h"
#include <string>
#include <fstream>
#include <vector>
using std::vector;

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        template<bool tTestHTheorem>
        class ELBM : public CollisionOperator
        {
          public:
            static void createAlphaArray(const size_t size);

            static void getSiteValues(const distribn_t* f,
                                      distribn_t &density,
                                      distribn_t &v_x,
                                      distribn_t &v_y,
                                      distribn_t &v_z,
                                      distribn_t* f_eq,
                                      const site_t index);

            static void getBoundarySiteValues(const distribn_t* f,
                                              const distribn_t &density,
                                              const distribn_t &v_x,
                                              const distribn_t &v_y,
                                              const distribn_t &v_z,
                                              distribn_t* f_eq,
                                              const site_t index);

            static void doPostCalculations(const distribn_t* f,
                                           const geometry::LatticeData* bLatDat,
                                           const site_t index);

            static distribn_t getOperatorElement(distribn_t &f_i,
                                                 distribn_t &f_eq_i,
                                                 const LbmParameters* iLbmParams);

            static void printHViolations();

          private:
            static double* alpha;
            static size_t currentAlphaIndex;
            static double HpreCollision;

            static vector<double> HViolationValues;
            static vector<unsigned int> HViolationTypes;

            static double getAlpha(const distribn_t* lFOld,
                                   const distribn_t* lFEq,
                                   double prevAlpha);
        };

        template<bool tTestHTheorem>
        double* ELBM<tTestHTheorem>::alpha;
        template<bool tTestHTheorem>
        size_t ELBM<tTestHTheorem>::currentAlphaIndex;
        template<bool tTestHTheorem>
        double ELBM<tTestHTheorem>::HpreCollision;

        template<bool tTestHTheorem>
        vector<double> ELBM<tTestHTheorem>::HViolationValues(0);
        template<bool tTestHTheorem>
        vector<unsigned int> ELBM<tTestHTheorem>::HViolationTypes(0);

        template<bool tTestHTheorem>
        void ELBM<tTestHTheorem>::createAlphaArray(const size_t size)
        {
          alpha = new double[size];
          for (size_t i = 0; i < size; i++)
            alpha[i] = 2.0;
        }

        template<bool tTestHTheorem>
        void ELBM<tTestHTheorem>::getSiteValues(const distribn_t* f,
                                                distribn_t &density,
                                                distribn_t &v_x,
                                                distribn_t &v_y,
                                                distribn_t &v_z,
                                                distribn_t* f_eq,
                                                const site_t index)
        {
          currentAlphaIndex = index;
          D3Q15::CalculateEntropicDensityVelocityFEq(f, density, v_x, v_y, v_z, f_eq);
          alpha[index] = getAlpha(f, f_eq, alpha[index]);

          if (tTestHTheorem)
          {
            HFunction HFunc(f, NULL);
            HpreCollision = HFunc.evaluate();
          }
        }

        template<bool tTestHTheorem>
        void ELBM<tTestHTheorem>::getBoundarySiteValues(const distribn_t* f,
                                                        const distribn_t &density,
                                                        const distribn_t &v_x,
                                                        const distribn_t &v_y,
                                                        const distribn_t &v_z,
                                                        distribn_t* f_eq,
                                                        const site_t index)
        {
          currentAlphaIndex = index;
          D3Q15::CalculateEntropicFeq(density, v_x, v_y, v_z, f_eq);
          alpha[index] = getAlpha(f, f_eq, alpha[index]);

          if (tTestHTheorem)
          {
            HFunction HFunc(f, NULL);
            HpreCollision = HFunc.evaluate();
          }
        }

        template<bool tTestHTheorem>
        void ELBM<tTestHTheorem>::doPostCalculations(const distribn_t* f,
                                                     const geometry::LatticeData* bLatDat,
                                                     const site_t index)
        {

          if (tTestHTheorem)
          {
            HFunction HFunc(f, NULL);
            double dH = HFunc.evaluate() - HpreCollision;

            // Cannot check change in H against zero as there is a possibility of round-off errors
            if (dH > 1.0E-10)
            {
              HViolationValues.push_back(dH);
              HViolationTypes.push_back(bLatDat->GetCollisionType(bLatDat->GetSiteData(index)));
            }
          }

        }

        template<bool tTestHTheorem>
        void ELBM<tTestHTheorem>::printHViolations()
        {
          if (tTestHTheorem)
          {
            if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
              return;

            std::string CollisionType;

            char filename[40];
            sprintf(filename,
                    "HViolations%d",
                    topology::NetworkTopology::Instance()->GetLocalRank());
            std::ofstream fout(filename);

            for (unsigned int i = 0; i < HViolationValues.size(); i++)
            {
              switch (HViolationTypes[i])
              {
                case FLUID:
                  CollisionType = "MidFluid";
                  break;
                case EDGE:
                  CollisionType = "Wall";
                  break;
                case INLET:
                  CollisionType = "Inlet";
                  break;
                case OUTLET:
                  CollisionType = "Outlet";
                  break;
                case (INLET | EDGE):
                  CollisionType = "InletWall";
                  break;
                case (OUTLET | EDGE):
                  CollisionType = "OutletWall";
                  break;
                default:
                  ;
              }

              fout << "!!!H theorem violated!!! dH = " << HViolationValues[i]
                  << "\tCollision type: " << CollisionType << std::endl;
            }

            fout.close();
          }
        }

        // Also updates lFEq_i to be lFNeq_i
        template<bool tTestHTheorem>
        distribn_t ELBM<tTestHTheorem>::getOperatorElement(distribn_t &f_i,
                                                           distribn_t &f_eq_i,
                                                           const LbmParameters* iLbmParams)
        {
          return (alpha[currentAlphaIndex] * iLbmParams->Beta * (f_eq_i = f_i - f_eq_i));
        }

        template<bool tTestHTheorem>
        double ELBM<tTestHTheorem>::getAlpha(const distribn_t* lF,
                                             const distribn_t* lFEq,
                                             double prevAlpha)
        {
          return 2.0;

          bool small = true;

          for (unsigned int i = 0; i < D3Q15::NUMVECTORS; i++)
          {
            // Papers suggest f_eq - f < 0.001 or (f_eq - f)/f < 0.01 for the point to have approx alpha = 2
            // Accuracy can change depending on stability requirements, because the more NR evaluations it skips
            // the more of the simulation is in the LBGK limit.
            if (fabs( (lFEq[i] - lF[i]) / lF[i]) > 1.0E-2)
            {
              small = false;
              break;
            }
          }

          HFunction HFunc(lF, lFEq);

          if (small)
          {
            double alphaLower = 1.8, HLower;
            double alphaHigher = 2.2, HHigher;

            HFunc(alphaLower, HLower);
            HFunc(alphaHigher, HHigher);

            // At the moment this decision is based on a few quick test cases
            // If the root is not near 2.0 then f - f_eq is too small to make the value of alpha matter
            // Chosen to return 2.0 as that is the LBGK case
            if (HLower * HHigher >= 0.0)
              return 2.0;

            return (hemelb::util::NumericalMethods::Brent(&HFunc, alphaLower, alphaHigher, 1.0E-3));
          }
          else
          {
            // This is in case previous Alpha was calculated to be zero (does happen occasionally if f_eq - f is small
            prevAlpha = (prevAlpha < 1.8
              ? 2.0
              : prevAlpha);

            // Accuracy is set to 1.0E-3 as this works for difftest.
            return (hemelb::util::NumericalMethods::NewtonRaphson(&HFunc, prevAlpha, 1.0E-3));
          }

        }

      }
    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ELBM_H */
