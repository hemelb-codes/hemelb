
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONDIRICHLET_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONDIRICHLET_H

#include "lb/streamers/AdvectionDiffusionStreamerTypeFactory.h"
#include "lb/streamers/AdvectionDiffusionDirichletDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"
#include "lb/streamers/GuoZhengShiDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<class CollisionType>
      struct AdvectionDiffusionDirichlet
      {
          typedef StentWallStreamerTypeFactory<CollisionType, AdvectionDiffusionDirichletDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionDirichletSBB
      {
          typedef StentWallVesselWallStreamerTypeFactory<CollisionType, AdvectionDiffusionDirichletDelegate<CollisionType> ,
              SimpleBounceBackDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionDirichletBFL
      {
          typedef StentWallVesselWallStreamerTypeFactory<CollisionType, AdvectionDiffusionDirichletDelegate<CollisionType> ,
              BouzidiFirdaousLallemandDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionDirichletGZS
      {
          typedef StentWallVesselWallStreamerTypeFactory<CollisionType, AdvectionDiffusionDirichletDelegate<CollisionType> ,
              GuoZhengShiDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionDirichletIoletSBB
      {
          typedef StentWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionDirichletDelegate<CollisionType> ,
              SimpleBounceBackDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionDirichletIoletBFL
      {
          typedef StentWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionDirichletDelegate<CollisionType> ,
              BouzidiFirdaousLallemandDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionDirichletIoletGZS
      {
          typedef StentWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionDirichletDelegate<CollisionType> ,
              GuoZhengShiDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionDirichletSBBIoletWall
      {
          typedef StentWallVesselWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionDirichletDelegate<CollisionType> ,
              SimpleBounceBackDelegate<CollisionType> , SimpleBounceBackDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionDirichletBFLIoletWall
      {
          typedef StentWallVesselWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionDirichletDelegate<CollisionType> ,
              BouzidiFirdaousLallemandDelegate<CollisionType> , BouzidiFirdaousLallemandDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionDirichletGZSIoletWall
      {
          typedef StentWallVesselWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionDirichletDelegate<CollisionType> ,
              GuoZhengShiDelegate<CollisionType> , GuoZhengShiDelegate<CollisionType> > Type;
      };


    }
  }
}

#endif /* HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONDIRICHLET_H */
