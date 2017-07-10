
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONBOUNDARYCONDITIONS_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONBOUNDARYCONDITIONS_H

#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/AdvectionDiffusionDirichletDelegate.h"
#include "lb/streamers/AdvectionDiffusionNeumannDelegate.h"
#include "lb/streamers/AdvectionDiffusionStreamerTypeFactory.h"

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
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONBOUNDARYCONDITIONS_H */