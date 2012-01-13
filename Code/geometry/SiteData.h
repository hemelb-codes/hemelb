#ifndef HEMELB_GEOMETRY_SITEDATA_H
#define HEMELB_GEOMETRY_SITEDATA_H

#include "units.h"

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
        SiteData(sitedata_t);
        SiteData(const SiteData& other);
        virtual ~SiteData();

        bool IsEdge() const;
        bool IsSolid() const;
        unsigned GetCollisionType() const;
        SiteType GetSiteType() const;
        int GetBoundaryId() const;
        bool HasBoundary(Direction direction) const;

        /**
         * This function should only be used for debugging.
         */
        sitedata_t GetRawValue() const;

        static const sitedata_t SITE_TYPE_BITS = 2U;
        static const sitedata_t BOUNDARY_CONFIG_BITS = 26U;
        static const sitedata_t BOUNDARY_DIR_BITS = 4U;
        static const sitedata_t BOUNDARY_ID_BITS = 10U;

        static const sitedata_t BOUNDARY_CONFIG_SHIFT = SITE_TYPE_BITS;
        static const sitedata_t BOUNDARY_DIR_SHIFT = BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS;
        static const sitedata_t BOUNDARY_ID_SHIFT = BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS;

        // Comments show the bit patterns.
        static const sitedata_t SITE_TYPE_MASK = ( ((sitedata_t) 1 << SITE_TYPE_BITS) - 1);
        // 0000 0000  0000 0000  0000 0000  0000 0011
        // These give the *_TYPE in geometry/LatticeData.h

        static const sitedata_t BOUNDARY_CONFIG_MASK = ( ((sitedata_t) 1 << BOUNDARY_CONFIG_BITS)
            - 1) << BOUNDARY_CONFIG_SHIFT;
        // 0000 0000  0000 0000  1111 1111  1111 1100
        // These bits are set if the lattice vector they correspond to takes one to a solid site
        // The following hex digits give the index into LatticeSite.neighbours
        // ---- ----  ---- ----  DCBA 9876  5432 10--

        static const sitedata_t BOUNDARY_DIR_MASK = ( ((sitedata_t) 1 << BOUNDARY_DIR_BITS) - 1)
            << BOUNDARY_DIR_SHIFT;
        // 0000 0000  0000 1111  0000 0000  0000 0000
        // No idea what these represent. As far as I can tell, they're unused.

        static const sitedata_t BOUNDARY_ID_MASK = ( ((sitedata_t) 1 << BOUNDARY_ID_BITS) - 1)
            << BOUNDARY_ID_SHIFT;
        // 0011 1111  1111 0000  0000 0000  0000 0000
        // These bits together give the index of the inlet/outlet/wall in the output XML file

        static const sitedata_t PRESSURE_EDGE_MASK = (sitedata_t) 1
            << (BOUNDARY_ID_BITS + BOUNDARY_ID_SHIFT + 1);
        // 1000 0000  0000 0000  0000 0000  0000 0000

      private:
        sitedata_t value;
    };
  }
}

#endif /* HEMELB_GEOMETRY_SITEDATA_H */
