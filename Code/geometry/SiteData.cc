// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "constants.h"
#include "geometry/SiteData.h"

namespace hemelb
{
  namespace geometry
  {
    /**
     * Constructor for a SiteData object
     *
     * NOTE: that this assumes a single fluid site will be next to at most one inlet or outlet.
     * This is a safe assumption in simple geometries but may not always be the case.
     *
     * @param readResult
     */
    SiteData::SiteData(const GeometrySite& readResult)
    {
      if (!readResult.isFluid)
      {
        data = SOLID_TYPE;
        boundaryIntersection = 0;
        ioletIntersection = 0;
      }
      else
      {
        int ioletId = 0;
        boundaryIntersection = 0;
        ioletIntersection = 0;

        bool hadInlet = false;
        bool hadOutlet = false;

        // Iterate over each direction
        for (Direction direction = 1; direction <= readResult.links.size(); ++direction)
        {
          // Get the link
          const GeometrySiteLink& link = readResult.links[direction - 1];

          // If it's a wall link, set the bit for this direction
          if (link.type == GeometrySiteLink::WALL_INTERSECTION)
          {
            boundaryIntersection |= 1 << (direction - 1);
          }
          else if (link.type == GeometrySiteLink::INLET_INTERSECTION || link.type
              == GeometrySiteLink::OUTLET_INTERSECTION)
          {
            ioletIntersection |= 1 << (direction - 1);
          }

          // If it's an inlet, take the IOlet id
          if (link.type == GeometrySiteLink::INLET_INTERSECTION)
          {
            ioletId = link.ioletId;
            hadInlet = true;
          }
          // Ditto if it's an outlet.
          else if (link.type == GeometrySiteLink::OUTLET_INTERSECTION)
          {
            ioletId = link.ioletId;
            hadOutlet = true;
          }
        }

        SiteType type = hadInlet
          ? INLET_TYPE
          : hadOutlet
            ? OUTLET_TYPE
            : FLUID_TYPE;

        data = (ioletId << BOUNDARY_ID_SHIFT) + type;
      }
    }

    SiteData::SiteData(const SiteData& other) :
      boundaryIntersection(other.boundaryIntersection), ioletIntersection(other.ioletIntersection), data(other.data)
    {
    }

    SiteData::SiteData() :
        boundaryIntersection(0), data(0)
    {
    }

    SiteData::~SiteData()
    {
    }

    bool SiteData::IsEdge() const
    {
      return boundaryIntersection != 0;
    }

    bool SiteData::IsSolid() const
    {
      return GetSiteType() == SOLID_TYPE;
    }

    unsigned SiteData::GetCollisionType() const
    {
      if (data == FLUID_TYPE && boundaryIntersection == 0)
      {
        return FLUID;
      }

      SiteType boundary_type = GetSiteType();

      if (boundary_type == FLUID_TYPE)
      {
        return EDGE;
      }
      if (boundaryIntersection == 0)
      {
        if (boundary_type == INLET_TYPE)
        {
          return INLET;
        }
        else
        {
          return OUTLET;
        }
      }
      else
      {
        if (boundary_type == INLET_TYPE)
        {
          return INLET | EDGE;
        }
        else
        {
          return OUTLET | EDGE;
        }
      }
    }

    SiteType SiteData::GetSiteType() const
    {
      return (SiteType) (data & SITE_TYPE_MASK);
    }

    int SiteData::GetBoundaryId() const
    {
      return (data & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
    }

    bool SiteData::HasBoundary(Direction direction) const
    {
      unsigned mask = 1U << (direction - 1);
      return (boundaryIntersection & mask) != 0;
    }

    bool SiteData::HasIolet(Direction direction) const
    {
      unsigned mask = 1U << (direction - 1);
      return (ioletIntersection & mask) != 0;
    }

    uint32_t SiteData::GetIoletIntersectionData() const
    {
      return ioletIntersection;
    }

    uint32_t SiteData::GetIntersectionData() const
    {
      return boundaryIntersection;
    }

    uint32_t SiteData::GetOtherRawData() const
    {
      return data;
    }

  }

}
