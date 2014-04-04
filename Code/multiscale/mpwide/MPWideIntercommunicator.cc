//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

/*
 * MPWideIntercommunicator.cc
 *
 *  Created on: 7 Nov 2012
 *      Author: Derek
 */

#include "multiscale/mpwide/MPWideIntercommunicator.h"
#include "net/IOCommunicator.h"
#include <MPWide.h>
#include <cstring>

namespace hemelb
{
  namespace multiscale
  {
    // Initialize the static initialized variable to false.
    bool MPWideIntercommunicator::mpwideInitialized = false;

    MPWideIntercommunicator::MPWideIntercommunicator(std::map<std::string, double> & buffer,
                                                     std::map<std::string, bool> &orchestration,
                                                     std::string configFilePathIn) :
        isCommsProc(net::IOCommunicator::Instance()->OnIORank()),
            configFilePath(configFilePathIn), recv_icand_data_size(0), send_icand_data_size(0),
            doubleContents(buffer), currentTime(0), orchestration(orchestration), channelCount(0)
    {
    }

    void MPWideIntercommunicator::Initialize()
    {
      if (isCommsProc)
      {
        log::Logger::Log<log::Info, log::Singleton>("Initializing MPWide.");

        // 1. Read the file with MPWide settings.
        std::vector < std::string > hosts;
        std::vector<int> server_side_ports;

        ReadInputFile(configFilePath.c_str(), hosts, server_side_ports);

        log::Logger::Log<log::Debug, log::Singleton>("MPWide input file read: base port is %i",
                                                     server_side_ports[0]);

        // 2. Initialize MPWide.
        MPW_Init(&hosts.front(), &server_side_ports.front(), channelCount);
      }
    }

    void MPWideIntercommunicator::ShareInitialConditions()
    {
      // If not yet initialized, initialize MPWide.
      if (!mpwideInitialized)
      {
        mpwideInitialized = true;
        Initialize();
      }

      if (isCommsProc)
      {
        // 1. Obtain and exchange shared data sizes.
        send_icand_data_size = GetRegisteredObjectsSize(registeredObjects);
        recv_icand_data_size = ExchangeICandDataSize(send_icand_data_size);

        hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("PRE-MALLOC, icand sizes are: %i (send) %i (recv)",
                                                                              send_icand_data_size,
                                                                              recv_icand_data_size);

        // 2. Allocate exchange buffers. We do this only once at initialization.
        ICandRecvDataPacked.resize(recv_icand_data_size);
        ICandSendDataPacked.resize(send_icand_data_size);
      }

      // Update the time and perform an initial exchange with the multiscale.
      doubleContents["shared_time"] = 0.0;
      ExchangeWithMultiscale();
    }

    /* This is run at the start of every time step in the main HemeLB simulation. */
    bool MPWideIntercommunicator::DoMultiscale(double new_time)
    {
      // 1. Update the shared time, if we should take a time step.
      bool shouldAdvance = ShouldAdvance();
      if (shouldAdvance)
      {
        UpdateSharedTime(new_time);
      }

      // 2. Exchange ICands with the other code.
      ExchangeWithMultiscale();

      // 3. Return the bool telling HemeLB whether to perform a timestep.
      return shouldAdvance;
    }

    void MPWideIntercommunicator::ExchangeWithMultiscale()
    {
      // 1. Pack/Serialize local shared data.
      hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Beginning exchange with multiscale");
      SerializeRegisteredObjects(&ICandSendDataPacked.front(), registeredObjects);

      // 2. Exchange serialized shared data.
      hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Exchanging packaged data");
      ExchangePackages(&ICandSendDataPacked.front(), &ICandRecvDataPacked.front());

      // 3. Unpack and merged the two serialized shared data copies.
      hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Unpacking and merging received data");
      UnpackReceivedData(registeredObjects, &ICandRecvDataPacked.front());

      hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Exchange with multiscale completed");
    }

    /* TODO: Only public for unit-testing. */
    void MPWideIntercommunicator::UnitTestIncrementSharedTime()
    {
      if (isCommsProc)
      {
        if (currentTime >= doubleContents["shared_time"])
        {
          doubleContents["shared_time"] += 1.0;
          printf("shared time set to %lf\n", doubleContents["shared_time"]);
        }
      }
    }

    void MPWideIntercommunicator::UpdateSharedTime(double new_time)
    {
      currentTime = new_time;
      doubleContents["shared_time"] = new_time;
    }

    bool MPWideIntercommunicator::ShouldAdvance()
    {
      return doubleContents["shared_time"] >= currentTime;
    }

    size_t MPWideIntercommunicator::GetTypeSize(RuntimeType type)
    {
      if (type == RuntimeTypeTraits::GetType<double>())
      {
        return sizeof(double);
      }
      if (type == RuntimeTypeTraits::GetType<int>())
      {
        return sizeof(int);
      }
      if (type == RuntimeTypeTraits::GetType<int64_t>())
      {
        return sizeof(int64_t);
      }

      hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::OnePerCore>("Error in GetTypeSize(): The RuntimeType is not recognized.");
      exit(-1);
      return -1;
    }

    FILE *MPWideIntercommunicator::checkingFileOpen(const char *path, const char *mode)
    {
      // Open the file
      FILE *fp = fopen(path, mode);

      // If there was a problem, print the error and exit
      if (fp == NULL)
      {
        perror(path);
        exit(EXIT_FAILURE);
      }

      // If there wasn't a problem, return the file handle.
      return fp;
    }

    void MPWideIntercommunicator::ReadInputFile(const char* socketsFilePath,
                                                std::vector<std::string>& url,
                                                std::vector<int>& serverSidePorts)
    {
      // Open the file
      FILE* socketsFile = checkingFileOpen(socketsFilePath, "r");

      /**
       * Read the number of channels and the number of hosts.
       */
      int numberOfHosts;
      int d = fscanf(socketsFile, "%d%d", &channelCount, &numberOfHosts);

      if (d > -1000)
      {
        log::Logger::Log<log::Warning, log::OnePerCore>("number of hosts: %i, number of channels: %i",
                                                        numberOfHosts,
                                                        channelCount);
      }

      /* Creating channels array */
      channels.resize(channelCount);
      for (int i = 0; i < channelCount; i++)
      {
        channels[i] = i;
      }

      // Track how many streams have so far been encountering (for setting port numbers)
      int streamsSeenSoFar = 0;

      // For each host
      for (int j = 0; j < numberOfHosts; j++)
      {
        // Read in the host name, how many streams it has, and its base port number.
        char host[256];
        int basePort;
        int numberOfStreamsOnHost;

        fscanf(socketsFile, "%s%d%d", host, &basePort, &numberOfStreamsOnHost);
        log::Logger::Log<log::Warning, log::OnePerCore>("host: %i, base channel: %i, num_streams: %i",
                                                        j,
                                                        streamsSeenSoFar,
                                                        numberOfStreamsOnHost);
        streamsSeenSoFar += numberOfStreamsOnHost;

        // For every stream from the host, add the host and port to the list of urls and ports.
        for (int i = 0; i < numberOfStreamsOnHost; i++)
        {
          url.push_back((std::string) host);
          int port = basePort + i;
          serverSidePorts.push_back(port);
        }
      }

      // Close the file.
      fclose(socketsFile);
    }

    /**
     *  Pack/Serialize local shared data
     *  TODO: include Endian conversion in the future
     **/
    void MPWideIntercommunicator::SerializeRegisteredObjects(char *sendDataPointer,
                                                             ContentsType registeredIcands)
    {
      // REMINDER: ContentsType = std::map<Intercommunicand *, std::pair<IntercommunicandTypeT *, std::string> >
      if (!isCommsProc)
      {
        return;
      }

      // Iterate over registered intercommunicands
      for (ContentsType::iterator icandProperties = registeredIcands.begin();
          icandProperties != registeredIcands.end(); icandProperties++)
      {
        // Dereference the iterator
        hemelb::multiscale::Intercommunicand &icandContained = *icandProperties->first;
        IntercommunicandTypeT &icandType = *icandProperties->second.first;
        std::string &icandLabel = icandProperties->second.second;

        hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Name of Icand = %s",
                                                                              icandLabel.c_str());

        // For every field on the current intercommunicand...
        for (unsigned int sharedFieldIndex = 0;
            sharedFieldIndex < icandContained.SharedValues().size(); sharedFieldIndex++)
        {
          // Get info about it.
          std::string &sharedValueLabel = icandType.Fields()[sharedFieldIndex].first;
          size_t SharedValueSize = GetTypeSize(icandType.Fields()[sharedFieldIndex].second);

          // Get the buffers
          void* sendingDataBuffer = (void *) sendDataPointer;
          void* localBufferOfDataToSend =
              (void *) & (*icandContained.SharedValues()[sharedFieldIndex]);

          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("memcpy in PackObj: size = %d",
                                                                                SharedValueSize);
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("sendingDataBuffer %f",
                                                                                * (static_cast<double*>(sendingDataBuffer)));
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("localBufferOfDataToSend %f",
                                                                                localBufferOfDataToSend);

          // Copy the local data into the buffer to be sent, and advance the pointer into the send buffer.
          memcpy(sendingDataBuffer, localBufferOfDataToSend, SharedValueSize);
          sendDataPointer += SharedValueSize;

          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Shared value: %s %i %f",
                                                                                sharedValueLabel.c_str(),
                                                                                sharedFieldIndex,
                                                                                * (static_cast<double*>(localBufferOfDataToSend)));

          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("done.");
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Shared value:");
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("%s",
                                                                                sharedValueLabel.c_str());
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("%i",
                                                                                sharedFieldIndex);
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("%f",
                                                                                * (static_cast<double*>(localBufferOfDataToSend)));
        }
      }
    }

    /* Exchange Serialized shared object packages between processes. */
    void MPWideIntercommunicator::ExchangePackages(char* ICandSendDataPacked,
                                                   char* ICandRecvDataPacked)
    {
      if (isCommsProc && (send_icand_data_size > 0 || recv_icand_data_size > 0))
      {
        // If this is the comms proc and we're expecting to send or receive some data, communicate it over MPWide
        MPW_SendRecv(ICandSendDataPacked,
                     (long long int) send_icand_data_size,
                     ICandRecvDataPacked,
                     (long long int) recv_icand_data_size,
                     &channels.front(),
                     channelCount);
      }
    }

    void MPWideIntercommunicator::UnpackReceivedData(ContentsType registeredIcands,
                                                     char *receivedDataPointer)
    {
      // Only the comms proc performs unpacking
      if (isCommsProc)
      {
        // Iterated over the intercommunicands
        for (ContentsType::iterator intercommunicandData = registeredIcands.begin();
            intercommunicandData != registeredIcands.end(); intercommunicandData++)
        {
          // Get the contents of the iterator.
          hemelb::multiscale::Intercommunicand &icandContained = *intercommunicandData->first;
          IntercommunicandTypeT &icandType = *intercommunicandData->second.first;

          // For every shared field...
          for (unsigned int sharedFieldIndex = 0;
              sharedFieldIndex < icandContained.SharedValues().size(); sharedFieldIndex++)
          {
            // Get the size of the current field
            size_t SharedValueSize = GetTypeSize(icandType.Fields()[sharedFieldIndex].second);

            // Get the buffers from the received data (at the current position) and the local variables to copy into
            void* receivedBuffer = (void *) receivedDataPointer;
            void* sharedValueBuffer = (void *) & (*icandContained.SharedValues()[sharedFieldIndex]);

            // Copy the data across, and update our position into the data.
            memcpy(sharedValueBuffer, receivedBuffer, SharedValueSize);
            receivedDataPointer += SharedValueSize;

            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Shared value [%d] = %f",
                                                                                  sharedFieldIndex,
                                                                                  static_cast<double*>(sharedValueBuffer)[sharedFieldIndex]);
          }
        }
      }
    }

    size_t MPWideIntercommunicator::GetRegisteredObjectsSize(ContentsType registeredObjects)
    {
      size_t size = 0;

      // Iterate over all registered objects
      for (ContentsType::iterator intercommunicandData = registeredObjects.begin();
          intercommunicandData != registeredObjects.end(); intercommunicandData++)
      {
        // Get the contents of the iterator
        hemelb::multiscale::Intercommunicand &sharedObject = *intercommunicandData->first;
        IntercommunicandTypeT &icandType = *intercommunicandData->second.first;

        // For every field that's shared,
        for (unsigned int sharedFieldIndex = 0;
            sharedFieldIndex < sharedObject.SharedValues().size(); sharedFieldIndex++)
        {
          // add the fields size to the total
          size += GetTypeSize(icandType.Fields()[sharedFieldIndex].second);
        }
      }

      return size;
    }

    int64_t MPWideIntercommunicator::ExchangeICandDataSize(int64_t send_icand_data_size)
    {
      int64_t rsize = 0;
      if (isCommsProc)
      {
        // Only the comms proc uses MPWide to exchange the sizes of data.
        // MPWide only supports char buffers, hence the casts here.
        MPW_SendRecv( ((char *) &send_icand_data_size),
                     sizeof(int64_t),
                     ((char *) &rsize),
                     sizeof(int64_t),
                     &channels.front(),
                     1);
      }

      return rsize;
    }
  }
}
