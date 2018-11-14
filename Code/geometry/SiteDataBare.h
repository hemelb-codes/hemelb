
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_SITEDATABARE_H
#define HEMELB_GEOMETRY_SITEDATABARE_H

#include "geometry/GeometrySite.h"
#include "geometry/SiteType.h"
#include "lb/cuda_helper.h"
#include "units.h"

namespace hemelb
{
  namespace geometry
  {
    class SiteData
    {
      public:
        SiteData(const GeometrySite& siteReadResult);
        SiteData(const SiteData& other);
        SiteData(); //default constructor allows one to use operator[] for maps
        ~SiteData();

        bool IsWall() const;
        bool IsSolid() const;
        unsigned GetCollisionType() const;

        CUDA_HOST_DEVICE SiteType GetSiteType() const
        {
          return type;
        }

        CUDA_HOST_DEVICE SiteType& GetSiteType()
        {
          return type;
        }

        CUDA_HOST_DEVICE int GetIoletId() const
        {
          return ioletId;
        }

        CUDA_HOST_DEVICE int& GetIoletId()
        {
          return ioletId;
        }

        CUDA_HOST_DEVICE bool HasWall(Direction direction) const
        {
          unsigned mask = 1U << (direction - 1);
          return ((wallIntersection & mask) != 0) && (direction > 0);
        }

        CUDA_HOST_DEVICE bool HasIolet(Direction direction) const
        {
          unsigned mask = 1U << (direction - 1);
          return ((ioletIntersection & mask) != 0) && (direction > 0);
        }

        /**
         * These functions return internal representations and should only be used for debugging.
         */
        uint32_t GetIoletIntersectionData() const
        {
          return ioletIntersection;
        }

        uint32_t &GetIoletIntersectionData()
        {
          return ioletIntersection;
        }

        uint32_t GetWallIntersectionData() const
        {
          return wallIntersection;
        }

        uint32_t &GetWallIntersectionData()
        {
          return wallIntersection;
        }

      protected:
        /**
         * This is a bit mask for whether a wall is hit by links in each direction.
         */
        uint32_t wallIntersection;

        /**
         * This is a bit mask for whether an iolet is hit by links in each direction.
         */
        uint32_t ioletIntersection;

        SiteType type;
        int ioletId;
    };
  }
}

#endif /* HEMELB_GEOMETRY_SITEDATABARE_H */
