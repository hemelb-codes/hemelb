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
#include "mpiInclude.h"

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
          return MpiDataTypeTraits<T>::GetMpiDataType();
        }
    };

    namespace mpwide
    {
      bool mpwide_initialized;
    }

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

     == MALLOC VS NEW ==
     The use of malloc in this class looks unusual because we use new and delete everywhere else,
     but in this instance it's necessary because some compilers might have different padding usage when
     'new' is used, leading to erroneous data transfers when the coupling is performed across platforms.


     This is a very dumb example of an intercommunicator. It stores communicated examples in a string-keyed buffer
     By sharing the same buffer between multiple intercommunicator interfaces, one can mock the behaviour of
     interprocess communication.
     */
    class MPWideIntercommunicator : public hemelb::multiscale::Intercommunicator<MPWideRuntimeType>
    {
      public:
        MPWideIntercommunicator(std::map<std::string, double> & buffer,
                                std::map<std::string, bool> &orchestration,
                                std::string configFilePathIn);
        void Initialize();

        /* This is run at the start of the HemeLB simulation, after Initialize(). */
        void ShareInitialConditions();
        /* This is run at the start of every time step in the main HemeLB simulation. */
        bool DoMultiscale(double new_time);

        /* TODO: Only public for unit-testing. */
        void UnitTestIncrementSharedTime();

      private:
        std::string configFilePath;

        // Data sizes and pointers of shared data recv and send buffers.
        int64_t recv_icand_data_size;
        int64_t send_icand_data_size;
        char *ICandRecvDataPacked;
        char *ICandSendDataPacked;

        void UpdateSharedTime(double new_time);

        bool ExchangeWithMultiscale();

        bool ShouldAdvance();
        long long int GetTypeSize(RuntimeType type);
        FILE *my_fopen(const char *path, const char *mode);
        int ReadInputHead(const char* sockets_file);
        void ReadInputFile(const char* sockets_file, std::string* url, int* server_side_ports, int* num_channels);

        /* Pack/Serialize local shared data (this may include Endian conversion in the future) */
        void PackRegisteredObjects(char *ICandSendDataPacked, ContentsType registeredIcands);

        /* Exchange Serialized shared object packages between processes. */
        void ExchangePackages(char* ICandSendDataPacked, char* ICandRecvDataPacked);

        void UnpackAndMergeRegisteredObjects(ContentsType registeredIcands, char *ICandRecvDataPacked);
        int64_t GetRegisteredObjectsSize(ContentsType registeredObjects);
        int64_t ExchangeICandDataSize(int64_t send_icand_data_size);

        bool isCommsProc;

        std::map<std::string, double> &doubleContents;
        double currentTime;
        std::map<std::string, bool> & orchestration;

        /* MPWide specific parameters */
        std::string *hosts;
        int* server_side_ports;
        int* channels;
        int num_channels;

        /* Shared values administration. */
        //int  numRegisteredObjects;
        //int* registeredObjectsIdList;
    }
    ;
  }
}

#endif // HEMELB_MULTISCALE_MPWIDE_MPWIDEINTERCOMMUNICATOR_H
