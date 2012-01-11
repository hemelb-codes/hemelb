#ifndef HEMELB_GEOMETRY_SITEDATA_H
#define HEMELB_GEOMETRY_SITEDATA_H

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
        SiteData(unsigned);
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
        unsigned GetRawValue() const;

        static const unsigned int SITE_TYPE_BITS = 2U;
        static const unsigned int BOUNDARY_CONFIG_BITS = 14U;
        static const unsigned int BOUNDARY_DIR_BITS = 4U;
        static const unsigned int BOUNDARY_ID_BITS = 10U;

        static const unsigned int BOUNDARY_CONFIG_SHIFT = 2U; // SITE_TYPE_BITS;
        static const unsigned int BOUNDARY_DIR_SHIFT = 16U; // BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS;
        static const unsigned int BOUNDARY_ID_SHIFT = 20U; // BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS;

        // Comments show the bit patterns.
        static const unsigned int SITE_TYPE_MASK = ( (1 << SITE_TYPE_BITS) - 1);
        // 0000 0000  0000 0000  0000 0000  0000 0011
        // These give the *_TYPE in geometry/LatticeData.h

        static const unsigned int BOUNDARY_CONFIG_MASK = ( (1 << BOUNDARY_CONFIG_BITS) - 1)
            << BOUNDARY_CONFIG_SHIFT;
        // 0000 0000  0000 0000  1111 1111  1111 1100
        // These bits are set if the lattice vector they correspond to takes one to a solid site
        // The following hex digits give the index into LatticeSite.neighbours
        // ---- ----  ---- ----  DCBA 9876  5432 10--

        static const unsigned int BOUNDARY_DIR_MASK = ( (1 << BOUNDARY_DIR_BITS) - 1)
            << BOUNDARY_DIR_SHIFT;
        // 0000 0000  0000 1111  0000 0000  0000 0000
        // No idea what these represent. As far as I can tell, they're unused.

        static const unsigned int BOUNDARY_ID_MASK = ( (1 << BOUNDARY_ID_BITS) - 1)
            << BOUNDARY_ID_SHIFT;
        // 0011 1111  1111 0000  0000 0000  0000 0000
        // These bits together give the index of the inlet/outlet/wall in the output XML file

        static const unsigned int PRESSURE_EDGE_MASK = 1U << 31U;
        // 1000 0000  0000 0000  0000 0000  0000 0000

      private:
        unsigned value;
    };
  }
}

#endif /* HEMELB_GEOMETRY_SITEDATA_H */
