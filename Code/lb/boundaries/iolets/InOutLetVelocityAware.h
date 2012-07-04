/*
 * InOutLetVelocityAware.h
 *
 *  Created on: 26 Jun 2012
 *      Author: derek
 */

#ifndef HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETVELOCITYAWARE_H
#define HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETVELOCITYAWARE_H

#include "lb/boundaries/iolets/InOutLetMultiscale.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "extraction/IterableDataSource.h"
#include "log/Logger.h"
namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {

        /***
         * An inlet/outlet whose density is obtained from the multiscale intercommunicator
         * We envisage communication of velocity information outwards to other processes
         * the velocity SharedValue is a place-holder for this.
         * We set 0.1 as an arbitrary value for the output velocity for now.
         * We do not yet have an understanding of the necessary physics, or the necessary computational infrastructure,
         * for building boundary conditions which have access to their nearby velocity field.
         * The min and max pressure SharedValues are placeholders for information needed by the steering and visualisation code.
         * The steering and visualisation code requires minimum and maximum pressure values.
         */
        class InOutLetVelocityAware : public InOutLetMultiscale
        {
          public:
            InOutLetVelocityAware() :
                InOutLetMultiscale(), NDM(NULL)

            {
              /* Single argument constructor... is this an issue? */
            }
            /***
             * The shared values are registered through the initialiser-list syntactic sugar.
             */
            InOutLetVelocityAware(const InOutLetVelocityAware &other) :
                InOutLetMultiscale(other), NDM(other.NDM)
            {
              /* Add a velocity aware exchange here if needed?*/
            }

            void InitialiseNeighbouringSites(hemelb::geometry::neighbouring::NeighbouringDataManager *manager,
                                             std::map<unsigned int, site_t> invertedBoundaryList)
            {
              NDM = manager;
            }

            virtual ~InOutLetVelocityAware()
            {
            }

            //virtual void DoIO(TiXmlElement *parent, bool isLoading, configuration::SimConfig* simConfig);

            virtual InOutLet* Clone() const
            {
              InOutLetVelocityAware* copy = new InOutLetVelocityAware(*this);
              return copy;
            }
          protected:
            hemelb::geometry::neighbouring::NeighbouringDataManager *NDM;
        };
      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETVELOCITYAWARE_H */
