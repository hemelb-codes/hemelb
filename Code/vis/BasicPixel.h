// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_VIS_BASICPIXEL_H
#define HEMELB_VIS_BASICPIXEL_H

#include "net/mpi.h"

namespace hemelb
{
  namespace vis
  {
    class BasicPixel
    {
      public:
        BasicPixel();
        BasicPixel(int iIn, int jIn);

        int GetI() const;
        int GetJ() const;

        bool operator <(const BasicPixel& right) const;
        bool operator ==(const BasicPixel& right) const;

        /**
         * NB This is not virtual, but it will be shadowed by derived classes.
         * In the base case, we don't combine the pixels as they aren't associated
         * with any data.
         */
        void Combine(const BasicPixel& other) const
        {

        }

        /**
         * Produces an MPI Datatype object but doesn't commit it or manage its memory.
         * @return
         */
        static MPI_Datatype GetMPIType();

      protected:
        int i, j;
    };
  }

  namespace net
  {
    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::vis::BasicPixel>::RegisterMpiDataType();
  }
}

#endif /* HEMELB_VIS_BASICPIXEL_H */
