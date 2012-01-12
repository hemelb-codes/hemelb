#include "constants.h"
#include "geometry/SiteData.h"

namespace hemelb
{
  namespace geometry
  {
    SiteData::SiteData(sitedata_t value) :
      value(value)
    {

    }

    SiteData::SiteData(const SiteData& other) :
      value(other.value)
    {
    }

    SiteData::~SiteData()
    {

    }

    bool SiteData::IsEdge() const
    {
      return (value & PRESSURE_EDGE_MASK) != 0;
    }

    bool SiteData::IsSolid() const
    {
      return GetSiteType() == SOLID_TYPE;
    }

    unsigned SiteData::GetCollisionType() const
    {
      if (value == FLUID_TYPE)
      {
        return FLUID;
      }

      SiteType boundary_type = GetSiteType();

      if (boundary_type == FLUID_TYPE)
      {
        return EDGE;
      }
      if (! (value & PRESSURE_EDGE_MASK))
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
      return (SiteType) (value & SITE_TYPE_MASK);
    }

    int SiteData::GetBoundaryId() const
    {
      return (value & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
    }

    bool SiteData::HasBoundary(Direction direction) const
    {
      unsigned mask = 1U << (direction - 1);
      unsigned shiftedMask = mask << BOUNDARY_CONFIG_SHIFT;
      return (value & shiftedMask) != 0;
    }

    sitedata_t SiteData::GetRawValue() const
    {
      return value;
    }
  }
}
