// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_MPICOMMUNICATOR_H
#define HEMELB_NET_MPICOMMUNICATOR_H

//#include "units.h"
//#include "net/mpi.h"
#include "net/MpiError.h"
#include <boost/shared_ptr.hpp>

namespace hemelb
{
  namespace net
  {
    class MpiGroup;

    class MpiCommunicator
    {
      public:
        static MpiCommunicator World();

        /**
         * Constructor for an uninitialised communicator, equivalent to
         * MPI_COMM_NULL
         * @param communicator
         */
        MpiCommunicator();

        /**
         * Class has virtual methods so should have virtual d'tor.
         */
        virtual ~MpiCommunicator()
        {

        }
        /**
         * Assignment operator.
         * @param rhs
         * @return
         */
        //MpiCommunicator& operator=(const MpiCommunicator& rhs);

        /**
         * Returns the local rank on the communicator
         * @return
         */
        virtual int Rank() const;

        /**
         * Returns the size of the communicator (i.e. total number of procs involved).
         * @return
         */
        virtual int Size() const;

        /**
         * Creates a new communicator
         * @param Group which is a subset of the group of this communicator.
         * @return New communicator.
         */
        MpiCommunicator Create(const MpiGroup& grp) const;

        /**
         * Allow implicit casts to MPI_Comm
         * @return The underlying MPI communicator.
         */
        operator MPI_Comm() const
        {
          return *commPtr;
        }
        /**
         * Is this communicator valid? I.e. not equal to MPI_COMM_NULL.
         */
        operator bool() const
        {
          return (bool)commPtr;
        }
        /**
         * Returns the MPI group being used.
         * @return
         */
        MpiGroup Group() const;

        /**
         * Abort - should try to bring down all tasks, but no guarantees
         * @param errCode
         */
        void Abort(int errCode) const;

      private:
        /**
         * Constructor to get data needed from an MPI communicator
         * @param communicator
         */
        MpiCommunicator(MPI_Comm communicator, bool willOwn);

        boost::shared_ptr<MPI_Comm> commPtr;
    };

    bool operator==(const MpiCommunicator& comm1, const MpiCommunicator& comm2);
    bool operator!=(const MpiCommunicator& comm1, const MpiCommunicator& comm2);

  }
}

#endif /* HEMELB_NET_MPICOMMUNICATOR_H */
