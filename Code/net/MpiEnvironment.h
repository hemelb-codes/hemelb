
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MPIENVIRONMENT_H
#define HEMELB_NET_MPIENVIRONMENT_H

namespace hemelb
{
  namespace net
  {
    /**
     * Manage the MPI environment and provide query/abort functions.
     *
     * The first-constructed instance will be responsible for calling
     * MPI_Init (on construction) and MPI_Finalize (on destruction).
     *
     */
    class MpiEnvironment
    {
      public:
        /**
         * Initialize the MPI environment, if it has not already been
         * If it has, this does nothing.
         *
         * @param argc The number of arguments provided
         * @param argv The array of argument strings
         */
        MpiEnvironment(int& argc, char**& argv);
        /**
         * If this instance created the MPI env, shut it down.
         * Otherwise, does nothing.
         */
        ~MpiEnvironment();

        /**
         * Query if MPI is initialised
         * @return
         */
        static bool Initialized();
        /**
         * Query if MPI is finalised
         * @return
         */
        static bool Finalized();
        /**
         * Abort MPI. Exact behaviour is implementation defined.
         * Should not return.
         *
         * @param errorCode The error code to return to the environment.
         */
        static void Abort(int errorCode);

      private:
        // Make copy c'tor and copy assign private to ensure uncopyable.
        MpiEnvironment( const MpiEnvironment& );
        const MpiEnvironment& operator=( const MpiEnvironment& );

        /// Whether this instance called MPI_Init
        bool doesOwnMpi;
    };
  }
}
#endif //HEMELB_NET_MPIENVIRONMENT_H
