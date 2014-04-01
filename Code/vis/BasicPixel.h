// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
}

#endif /* HEMELB_VIS_BASICPIXEL_H */
