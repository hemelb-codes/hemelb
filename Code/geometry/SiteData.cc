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
#include "mpiInclude.h"
namespace hemelb
{
  template<>
  MPI_Datatype MpiDataTypeTraits<geometry::SiteType>::RegisterMpiDataType()
  {
    return MpiDataTypeTraits<int>::RegisterMpiDataType();
  }

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
        type = SOLID_TYPE;
        wallIntersection = 0;
        ioletIntersection = 0;
        ioletId = -1;
      }
      else
      {
        ioletId = -1;
        wallIntersection = 0;
        ioletIntersection = 0;

        bool hadInlet = false;
        bool hadOutlet = false;

        // Iterate over each (non-zero) direction
        for (Direction direction = 1; direction <= readResult.links.size(); ++direction)
        {
          // Get the link
          const GeometrySiteLink& link = readResult.links[direction - 1];

          // If it's a wall link, set the bit for this direction
          if (link.type == GeometrySiteLink::WALL_INTERSECTION)
          {
            wallIntersection |= 1 << (direction - 1);
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

        type = hadInlet
          ? INLET_TYPE
          : (hadOutlet
            ? OUTLET_TYPE
            : FLUID_TYPE);

      }
    }

    SiteData::SiteData(const SiteData& other) :
      wallIntersection(other.wallIntersection), ioletIntersection(other.ioletIntersection),
          type(other.type), ioletId(other.ioletId)
    {
    }

    SiteData::SiteData() :
      wallIntersection(0), ioletIntersection(0), type(SOLID_TYPE), ioletId(-1)
    {
    }

    SiteData::~SiteData()
    {
    }

    bool SiteData::IsWall() const
    {
      return wallIntersection != 0;
    }

    bool SiteData::IsSolid() const
    {
      return GetSiteType() == SOLID_TYPE;
    }

    unsigned SiteData::GetCollisionType() const
    {
      if (wallIntersection == 0)
      {
        // No solid wall intersections
        switch (type)
        {
          case FLUID_TYPE:
            return FLUID;

          case INLET_TYPE:
            return INLET;

          case OUTLET_TYPE:
            return OUTLET;
        }
      }
      else
      {
        // There are solid wall intersections
        switch (type)
        {
          case FLUID_TYPE:
            return WALL;

          case INLET_TYPE:
            return INLET | WALL;

          case OUTLET_TYPE:
            return OUTLET | WALL;
        }

      }
    }

    bool SiteData::HasWall(Direction direction) const
    {
      unsigned mask = 1U << (direction - 1);
      return (wallIntersection & mask) != 0;
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

    uint32_t SiteData::GetWallIntersectionData() const
    {
      return wallIntersection;
    }

  }

}
