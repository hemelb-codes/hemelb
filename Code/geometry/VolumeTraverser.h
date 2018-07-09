
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_VOLUMETRAVERSER_H
#define HEMELB_GEOMETRY_VOLUMETRAVERSER_H

#include "util/Vector3D.h"

namespace hemelb
{
  namespace geometry
  {
    /**
     * VolumeTraverser is used to sequentially traverse a 
     * 3D structure maintaining the index and Location 
     * within volume
     */
    class VolumeTraverser
    {
      public:
        virtual ~VolumeTraverser();

        util::Vector3D<site_t> GetCurrentLocation();

        void SetCurrentLocation(const util::Vector3D<site_t>& iLocation);

        site_t GetX();

        site_t GetY();

        site_t GetZ();

        site_t GetCurrentIndex() const;

        site_t GetIndexFromLocation(util::Vector3D<site_t> iLocation) const;

        //Increments the index by one and update the location accordingly
        //Returns true if successful or false if the whole volume has been
        //traversed
        bool TraverseOne();

        void IncrementX();
        void IncrementY();
        void IncrementZ();

        void DecrementX();
        void DecrementY();
        void DecrementZ();

        bool CurrentLocationValid();

        //Virtual methods which must be defined for correct traversal
        virtual site_t GetXCount() const = 0;
        virtual site_t GetYCount() const = 0;
        virtual site_t GetZCount() const = 0;

      protected:
        VolumeTraverser();

      private:
        util::Vector3D<site_t> mCurrentLocation;
        site_t mCurrentNumber;
    };
  }
}

#endif // HEMELB_GEOMETRY_VOLUMETRAVERSER_H
