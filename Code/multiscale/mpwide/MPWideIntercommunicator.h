// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_MULTISCALE_MPWIDE_MPWIDEINTERCOMMUNICATOR_H
#define HEMELB_MULTISCALE_MPWIDE_MPWIDEINTERCOMMUNICATOR_H

#include <unistd.h>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <vector>
#include <sstream>

#include "multiscale/Intercommunicator.h"
#include "net/mpi.h"

namespace hemelb
{
  namespace multiscale
  {
    /***
     * Example type traits structure, using the HemeLB implementation of MPI_TYPE traits.
     */
    struct MPWideRuntimeType
    {
        typedef MPI_Datatype RuntimeType;
        template<class T> static RuntimeType GetType()
        {
          return net::MpiDataTypeTraits<T>::GetMpiDataType();
        }
    };

    /**
     MPWide Intercommunicator class.
     This class provides an abstraction for the MPWide Communication Library in HemeLB.
     http://castle.strw.leidenuniv.nl/software/MPWide.html

     Dictionary:
     ICand = Intercommunicand

     Architectural Assumptions in this Intercommunicator:
     - We assume that shared value collections are constant in size throughout the simulation.
     + This will need to be changed if we introduce adaptive elements in the boundaries such
     as flexible vessel walls.

     - We assume that both HemeLB end points use the same type of ICand. These ICands can be
     different in size however if they contain vectors individually.
     + We currently think this limitation actually encourages writing proper ICands.

     This is a very dumb example of an intercommunicator. It stores communicated examples in a string-keyed buffer
     By sharing the same buffer between multiple intercommunicator interfaces, one can mock the behaviour of
     interprocess communication.
     */
    class MPWideIntercommunicator : public hemelb::multiscale::Intercommunicator<MPWideRuntimeType>
    {
      public:
        MPWideIntercommunicator(bool isCommsRank, std::map<std::string, double> & buffer,
                                std::map<std::string, bool> &orchestration,
                                std::string configFilePathIn);
        /** This is run at the start of the HemeLB simulation. */
        void ShareInitialConditions();
        /** This is run at the start of every time step in the main HemeLB simulation. */
        bool DoMultiscale(double new_time);

        /** TODO: Only public for unit-testing. */
        void UnitTestIncrementSharedTime();

      private:

        /**
         *  Initialises MPWide.
         */
        void Initialize();

        /**
         * Updates the shared value of time to a new value.
         * @param new_time
         */
        void UpdateSharedTime(double new_time);

        /**
         * Performs the multiscale exchange over MPWide.
         */
        void ExchangeWithMultiscale();

        /**
         * True if we should advance the current time.
         * @return
         */
        bool ShouldAdvance();

        /**
         * Get the size of the C++ type associated.
         * @param type
         * @return
         */
        size_t GetTypeSize(RuntimeType type);

        /**
         * Open file, checking for problems.
         * @param path
         * @param mode
         * @return
         */
        FILE *checkingFileOpen(const char *path, const char *mode);

        /**
         * Read the entirety of the MPWide config file
         * @param sockets_file
         * @param url
         * @param server_side_ports
         */
        void ReadInputFile(const char* sockets_file, std::vector<std::string>& url,
                           std::vector<int>& server_side_ports);

        /**
         *  Serialize local shared data (this may include Endian conversion in the future)
         **/
        void SerializeRegisteredObjects(char *ICandSendDataPacked, ContentsType registeredIcands);

        /**
         *  Exchange Serialized shared object packages between processes.
         **/
        void ExchangePackages(char* ICandSendDataPacked, char* ICandRecvDataPacked);

        /**
         * Unpack the received data from the intercommunicand
         * @param registeredIcands
         * @param ICandRecvDataPacked
         */
        void UnpackReceivedData(ContentsType registeredIcands, char *ICandRecvDataPacked);

        /**
         * Get the total size of the objects registered with the intercommunicand
         * @param registeredObjects
         * @return
         */
        size_t GetRegisteredObjectsSize(ContentsType registeredObjects);

        /**
         * Exchange data size of intercommunicand
         *
         * @param send_icand_data_size
         * @return Received data size
         */
        int64_t ExchangeICandDataSize(int64_t send_icand_data_size);

        /**
         * True if this is the proc that does all the communication over MPWide.
         */
        bool isCommsProc;

        /**
         *  The path to the MPWide config file.
         */
        std::string configFilePath;

        /**
         * Data sizes and pointers of shared data recv and send buffers.
         */
        int64_t recv_icand_data_size;
        int64_t send_icand_data_size;
        std::vector<char> ICandRecvDataPacked;
        std::vector<char> ICandSendDataPacked;

        /**
         * Map of string to double for the shared time over the intercommunicand
         */
        std::map<std::string, double> &doubleContents;

        /**
         * Current time.
         */
        double currentTime;

        /**
         * Map telling what things should be orchestrated.
         */
        std::map<std::string, bool> & orchestration;

        /**
         * Number of MPIWide channels.
         **/
        int channelCount;
        /**
         * The ids of the MPWide channels.
         */
        std::vector<int> channels;

        /**
         * True if MPWide has been initialised.
         */
        static bool mpwideInitialized;
    };
  }
}

#endif // HEMELB_MULTISCALE_MPWIDE_MPWIDEINTERCOMMUNICATOR_H
