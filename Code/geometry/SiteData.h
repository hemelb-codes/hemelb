#ifndef HEMELB_GEOMETRY_SITEDATA_H
#define HEMELB_GEOMETRY_SITEDATA_H

#include "units.h"
#include "geometry/GeometrySite.h"

namespace hemelb
{
  namespace geometry
  {
    enum SiteType
    {
      // These must be consistent with the setup tool
      SOLID_TYPE = 0U,
      FLUID_TYPE = 1U,
      INLET_TYPE = 2U,
      OUTLET_TYPE = 3U
    };

    class SiteData
    {
      public:
        SiteData(const GeometrySite& siteReadResult);
        SiteData(const SiteData& other);
        SiteData() :
          boundaryIntersection(0), ioletIntersection(0), data(0)
        {
        }//default constructor allows one to use operator[] for maps
        virtual ~SiteData();

        bool IsEdge() const;
        bool IsSolid() const;
        unsigned GetCollisionType() const;
        SiteType GetSiteType() const;
        int GetBoundaryId() const;
        bool HasBoundary(Direction direction) const;
        bool HasIolet(Direction direction) const;

        /**
         * These functions return internal representations and should only be used for debugging.
         */
        uint32_t GetIntersectionData() const;
        uint32_t &GetIntersectionData()
        {
          return boundaryIntersection;
        }
        uint32_t GetOtherRawData() const;
        uint32_t &GetOtherRawData()
        {
          return data;
        }
        static const uint32_t SITE_TYPE_BITS = 2U;
        /**
         * Arbitrarily set this to 20 bits.
         */
        static const uint32_t BOUNDARY_ID_BITS = 20U;

        static const uint32_t SITE_TYPE_MASK = (1 << SITE_TYPE_BITS) - 1;

        static const sitedata_t BOUNDARY_ID_MASK = ( (1 << (SITE_TYPE_BITS + BOUNDARY_ID_BITS)) - 1) - SITE_TYPE_MASK;
        static const sitedata_t BOUNDARY_ID_SHIFT = SITE_TYPE_BITS;

      protected:
        /**
         * This is a bit mask for whether a wall is hit by links in each direction.
         */
        uint32_t boundaryIntersection;

        /**
         * This is a bit mask for whether an iolet is hit by links in each direction.
         */
        uint32_t ioletIntersection;

        /**
         * This is a bit field that stores all other data associated with the site:
         * 2 bits for the site type.
         * The remaining 30 bits are for specifying the id of the boundary.
         */
        uint32_t data;
    };
  }
}

#endif /* HEMELB_GEOMETRY_SITEDATA_H */
