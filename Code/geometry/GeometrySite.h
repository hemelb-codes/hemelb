#ifndef HEMELB_GEOMETRY_GEOMETRYSITE_H
#define HEMELB_GEOMETRY_GEOMETRYSITE_H

#include <vector>
#include "constants.h"
#include "units.h"
#include "geometry/GeometrySiteLink.h"

namespace hemelb
{
  namespace geometry
  {
    /***
     * Model of the data for a site, as contained within a geometry file.
     * this data will be broken up and placed in various arrays in hemelb::Geometry::LatticeData
     *
     * Note that this should be able to be returned by copy (as we sometimes do) so be careful about
     * using heap-allocated data in this struct.
     */
    struct GeometrySite
    {
      public:
        //! Basic constructor for solid and fluid sites.
        GeometrySite(bool siteIsFluid) :
            targetProcessor(siteIsFluid ?
              -1 :
              BIG_NUMBER2), isFluid(siteIsFluid)
        {
        }

        GeometrySite(const GeometrySite& other) :
            targetProcessor(other.targetProcessor), isFluid(other.isFluid), links(other.links)
        {

        }

        //! Processor on which to perform lattice-Boltzmann for the site.
        proc_t targetProcessor;

        //! True iff the site is fluid, i.e. it is within the geometry and we will be doing
        //! lattice-Boltzmann with it.
        bool isFluid;

        //! A vector of the link data for each direction in the lattice currently being used
        //! (NOT necessarily the same as the lattice used by the geometry file).
        std::vector<GeometrySiteLink> links;
    };
  }
}

#endif /* HEMELB_GEOMETRY_SITEREADRESULT_H */
