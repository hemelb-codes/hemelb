
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONBOUNDARYCONDITIONS_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONBOUNDARYCONDITIONS_H

#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/AdvectionDiffusionDirichletDelegate.h"
#include "lb/streamers/AdvectionDiffusionInletDirichletDelegate.h"
#include "lb/streamers/AdvectionDiffusionNeumannDelegate.h"
#include "lb/streamers/AdvectionDiffusionOutflowDelegate.h"
#include "lb/streamers/AdvectionDiffusionOutflowBounceBackDelegate.h"
#include "lb/streamers/AdvectionDiffusionStreamerTypeFactory.h"
#include "lb/streamers/AdvectionDiffusionOutflowStreamerTypeFactory.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLDirichlet
      {
          typedef AdvectionDiffusionWallStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLNeumann
      {
          typedef AdvectionDiffusionWallStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLIolet
      {
          typedef AdvectionDiffusionIoletStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLWallBFLIoletDirichlet
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLWallBFLIoletNeumann
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLWallSBBIoletDirichlet
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLWallSBBIoletNeumann
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBDirichlet
      {
          typedef AdvectionDiffusionWallStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBNeumann
      {
          typedef AdvectionDiffusionWallStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBIolet
      {
          typedef AdvectionDiffusionIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBWallBFLIoletDirichlet
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBWallBFLIoletNeumann
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBWallSBBIoletDirichlet
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBWallSBBIoletNeumann
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionOutflowIolet
      {
          typedef AdvectionDiffusionOutflowIoletStreamerTypeFactory<CollisionImpl, AdvectionDiffusionOutflowDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBWallOutflowIoletDirichlet
      {
          typedef AdvectionDiffusionWallOutflowIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionOutflowDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBWallOutflowIoletNeumann
      {
          typedef AdvectionDiffusionWallOutflowIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionOutflowDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLWallOutflowIoletDirichlet
      {
          typedef AdvectionDiffusionWallOutflowIoletStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionOutflowDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLWallOutflowIoletNeumann
      {
          typedef AdvectionDiffusionWallOutflowIoletStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionOutflowDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionInletDirichlet
      {
          typedef AdvectionDiffusionIoletStreamerTypeFactory<CollisionImpl, AdvectionDiffusionInletDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBWallInletDirichletDirichlet
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionInletDirichletDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBWallInletDirichletNeumann
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionInletDirichletDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLWallInletDirichletDirichlet
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionInletDirichletDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLWallInletDirichletNeumann
      {
          typedef AdvectionDiffusionWallIoletStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionInletDirichletDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionOutflowBounceBackIolet
      {
          typedef AdvectionDiffusionOutflowIoletStreamerTypeFactory<CollisionImpl, AdvectionDiffusionOutflowBounceBackDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBWallOutflowBounceBackIoletDirichlet
      {
          typedef AdvectionDiffusionWallOutflowIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionOutflowBounceBackDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionSBBWallOutflowBounceBackIoletNeumann
      {
          typedef AdvectionDiffusionWallOutflowIoletStreamerTypeFactory<CollisionImpl, SimpleBounceBackDelegate<CollisionImpl>, AdvectionDiffusionOutflowBounceBackDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLWallOutflowBounceBackIoletDirichlet
      {
          typedef AdvectionDiffusionWallOutflowIoletStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionOutflowBounceBackDelegate<CollisionImpl>, AdvectionDiffusionDirichletDelegate<CollisionImpl> > Type;
      };

      template<typename CollisionImpl>
      struct AdvectionDiffusionBFLWallOutflowBounceBackIoletNeumann
      {
          typedef AdvectionDiffusionWallOutflowIoletStreamerTypeFactory<CollisionImpl, BouzidiFirdaousLallemandDelegate<CollisionImpl>, AdvectionDiffusionOutflowBounceBackDelegate<CollisionImpl>, AdvectionDiffusionNeumannDelegate<CollisionImpl> > Type;
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONBOUNDARYCONDITIONS_H */
