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
#include "multiscale/Intercommunicator.h"
#include "mpiInclude.h"
#include "tinyxml.h"

/**
 MPWide Intercommunicator class.
 This class provides an abstraction for the MPWide Communication Library in HemeLB.
 http://castle.strw.leidenuniv.nl/software/MPWide.html
 **/

#include <algorithm>
#include <functional>
#include <vector>
#include <sstream>
/* TODO: make a good separation of test and production includes here. */
//#include "unittests/multiscale/MockMPWide.h" /* This is temporary! */
//#include "MPWide.h"
#include <unistd.h>
#include <cstdio>

/*
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

 */

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

    /***
     * This is a very dumb example of an intercommunicator
     * It stores communicated examples in a string-keyed buffer
     * By sharing the same buffer between multiple intercommunicator interfaces, one can mock the behaviour of interprocess communication.
     */
    namespace mpwide
    {
      bool mpwide_initialized = false;
      std::string mpwide_config_file = "../../../config_files/MPWSettings.cfg";
      bool mpwide_comm_proc = false;
    }

    class MPWideIntercommunicator : public hemelb::multiscale::Intercommunicator<MPWideRuntimeType>
    {
      public:
        MPWideIntercommunicator(std::map<std::string, double> & buffer, std::map<std::string, bool> &orchestration) :
            doubleContents(buffer), currentTime(0), orchestration(orchestration)
        {
          /* if(!hemelb::multiscale::mpwide::mpwide_initialized) {
           hemelb::multiscale::mpwide::mpwide_initialized = true;
           Initialize();
           }*/
        }

        void Initialize()
        {
          //TODO: This is a temporary hard-code.
          //The MPWide config file should be read from the HemeLB XML config file!

          if (topology::NetworkTopology::Instance()->GetLocalRank() == 0)
          {
            hemelb::multiscale::mpwide::mpwide_comm_proc = true;
          }

          if (hemelb::multiscale::mpwide::mpwide_comm_proc)
          {
            fprintf(stdout, "Running MPWide Initialize().\n");

            char cwd[1024];
            if (getcwd(cwd, sizeof (cwd)) != NULL)
            {
              fprintf(stdout, "Current [head] working dir: %s\n", cwd);
            }
            else
            {
              perror("getcwd() error");
            }

            num_channels = ReadInputHead(hemelb::multiscale::mpwide::mpwide_config_file.c_str());

            /* Creating channels array */
            channels = (int *) malloc(num_channels * sizeof(int));
            for (int i = 0; i < num_channels; i++)
            {
              channels[i] = i;
            }

            hosts = new std::string[num_channels];
            server_side_ports = (int *) malloc(num_channels * sizeof(int));

            if (getcwd(cwd, sizeof (cwd)) != NULL)
            {
              fprintf(stdout, "Current working dir: %s\n", cwd);
            }
            else
            {
              perror("getcwd() error");
            }

            // 1. Read the file with MPWide settings.
            ReadInputFile(hemelb::multiscale::mpwide::mpwide_config_file.c_str(),
                          hosts,
                          server_side_ports,
                          &num_channels);

            std::cout << "MPWide input file read: base port = " << server_side_ports << std::endl;

            // 2. Initializa MPWide.
            MPW_Init(hosts, server_side_ports, num_channels);
          }
        }

        /* This is run at the start of the HemeLB simulation, after Initialize(). */
        void ShareInitialConditions()
        {
          if (!hemelb::multiscale::mpwide::mpwide_initialized)
          {
            hemelb::multiscale::mpwide::mpwide_initialized = true;
            Initialize();
          }

          if (hemelb::multiscale::mpwide::mpwide_comm_proc)
          {
            // 1. Obtain and exchange shared data sizes.
            send_icand_data_size = GetRegisteredObjectsSize(registeredObjects);
            recv_icand_data_size = ExchangeICandDataSize(send_icand_data_size);
            //recv_icand_data_size = send_icand_data_size; /* TODO: DeMock! */

            //TODO: Add an offset table for the ICand data.

            std::cout << "PRE-MALLOC, icand sizes are: " << send_icand_data_size << "/" << recv_icand_data_size
                << std::endl;

            // 2. Allocate exchange buffers. We do this once at initialization,
            //    so that if it goes wrong, the program will crash timely.
            ICandRecvDataPacked = (char *) malloc(recv_icand_data_size);
            ICandSendDataPacked = (char *) malloc(send_icand_data_size);
          }
            doubleContents["shared_time"] = 0.0;
            ExchangeWithMultiscale(); //is this correct???

        }

        /* This is run at the start of every time step in the main HemeLB simulation. */
        bool DoMultiscale(double new_time)
        {
          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Entering DM");
          // 1. check whether this HemeLB instance was supposed to take a time step.
          if (ShouldAdvance())
          {

            // 2. Update the shared time.
            UpdateSharedTime(new_time);
          }
          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Calling EWM");
          // 3. Exchange ICands with the other code.
          ExchangeWithMultiscale();

          // 4. Return the bool telling HemeLB whether to perform a timestep.
          return ShouldAdvance();
        }

        /* TODO: Only public for unit-testing. */
        void UnitTestIncrementSharedTime()
        {
          if (hemelb::multiscale::mpwide::mpwide_comm_proc)
          {
            if (currentTime >= doubleContents["shared_time"])
            {
              doubleContents["shared_time"] += 1.0;
              printf("shared time set to %lf\n", doubleContents["shared_time"]);
            }
          }
        }

      private:

        // Data sizes and pointers of shared data recv and send buffers.
        int64_t recv_icand_data_size;
        int64_t send_icand_data_size;
        char *ICandRecvDataPacked;
        char *ICandSendDataPacked;

        void UpdateSharedTime(double new_time)
        {
          currentTime = new_time;
          doubleContents["shared_time"] = new_time;
        }

        bool ShouldAdvance()
        {
          return doubleContents["shared_time"] >= currentTime;
        }

        long long int GetTypeSize(RuntimeType type, hemelb::multiscale::BaseSharedValue & value)
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

          std::cerr << "Error in GetTypeSize(): The RuntimeType is not recognized." << std::endl;
          exit(-1);
          return -1;
        }

        FILE *my_fopen(const char *path, const char *mode)
        {

          FILE *fp;
          fp = fopen(path, mode);

          if (fp == NULL)
          {
            perror(path);
            exit(EXIT_FAILURE);
          }

          return fp;

        }

        int ReadInputHead(const char* sockets_file)
        {
          int nstream;
          int nhost;

          FILE *fin = my_fopen(sockets_file, "r");

          /* in the sockets_file, nstream indicates the total number of streams. */
          int d = fscanf(fin, "%d%d", &nstream, &nhost);

          if (d > -1000)
          {
            std::cerr << "nhost: " << nhost << ", nstream_total: " << nstream << std::endl;
          }

          fclose(fin);

          return nstream;
        }

        void ReadInputFile(const char* sockets_file, std::string* url, int* server_side_ports, int* num_channels)
        {

          int nstream;
          int nhost;
          char host[256];

          FILE *fin = my_fopen(sockets_file, "r");

          /* in the sockets_file, nstream indicates the total number of streams. */
          int d = fscanf(fin, "%d%d", &nstream, &nhost);

          /* for n sites, every site has two neighbours, and will establish
           * links with these neighbours. */
          //string hosts[nstream];
          //int ports[nstream];
          //int cports[nstream];
          //server_side_ports = (int *) malloc(nstream*sizeof(int));
          num_channels[0] = nstream;

          int nstream_host[nhost];
          int base_host_channel[nhost];
          int base_port[nhost];

          int streamcount_tmp = 0;
          int offset = 0;

          for (int j = 0; j < nhost; j++)
          {
            d = fscanf(fin, "%s%d%d", host, & (base_port[j]), & (nstream_host[j]));
            base_host_channel[j] = streamcount_tmp;
            streamcount_tmp += nstream_host[j];
            std::cerr << "host: " << j << ", base channel: " << base_host_channel[j] << ", num_streams: "
                << nstream_host[j] << std::endl;

            for (int i = 0; i < nstream_host[j]; i++)
            {
              int ii = offset + i;
              //std::cerr << "ii defined: " << ii << std::endl;
              url[ii] = (std::string) host;
              //std::cerr << "hosts defined: " << url[ii] << std::endl;
              int p = base_port[j] + i;
              //std::cerr << "p defined. Nstream is " << nstream << std::endl;
              server_side_ports[ii] = p;
              //cports[ii] = ports[ii] + nstream_host[j];
              //std::cerr << url[ii] << "\t" << server_side_ports[ii] << std::endl;
            }
            offset += nstream_host[j];
          }

          fclose(fin);

          //std::cout << "MPWide Settings File has been processed... " << std::endl;
        }

        /* Pack/Serialize local shared data (this may include Endian conversion in the future) */
        void PackRegisteredObjects(char *ICandSendDataPacked, ContentsType registeredObjects)
        {
          /* REMINDER:
           ContentsType = std::map<Intercommunicand *, std::pair<IntercommunicandTypeT *, std::string> >
           */

          int64_t offset = 0;

          if (hemelb::multiscale::mpwide::mpwide_comm_proc)
          {
            for (ContentsType::iterator intercommunicandData = registeredObjects.begin();
                intercommunicandData != registeredObjects.end(); intercommunicandData++)
            {
              hemelb::multiscale::Intercommunicand &sharedObject = *intercommunicandData->first;
              //std::string &label = intercommunicandData->second.second;
              IntercommunicandTypeT &resolver = *intercommunicandData->second.first;

              for (unsigned int sharedFieldIndex = 0; sharedFieldIndex < sharedObject.Values().size();
                  sharedFieldIndex++)
              {
                int64_t size = GetTypeSize(resolver.Fields()[sharedFieldIndex].second,
                                           *sharedObject.Values()[sharedFieldIndex]);
                void *buf1 = (void *) & (ICandSendDataPacked[offset]);
                void *buf2 = (void *) & (*sharedObject.Values()[sharedFieldIndex]);
                std::cout << "memcpy in PackObj: size = " << size << std::endl;
                memcpy(buf1, buf2, size);
                std::cout << "done." << std::endl;
                offset += size;
              }
            }
          }

        }

        /* Exchange Serialized shared object packages between processes. */
        void ExchangePackages(char* ICandSendDataPacked, char* ICandRecvDataPacked)
        {
          if (hemelb::multiscale::mpwide::mpwide_comm_proc && (send_icand_data_size > 0 || recv_icand_data_size > 0))
          {
            std::cout << "EXCHANGE PACKAGES, icand sizes are: " << send_icand_data_size << "/" << recv_icand_data_size
                << std::endl;
            MPW_SendRecv(ICandSendDataPacked,
                         (long long int) send_icand_data_size,
                         ICandRecvDataPacked,
                         (long long int) recv_icand_data_size,
                         channels,
                         num_channels);

            std::cout << "Value received (CHEAP DEBUG!): " << ((double*) ICandRecvDataPacked)[0] << std::endl;
          }
        }

        void UnpackAndMergeRegisteredObjects(ContentsType registeredObjects, char *ICandRecvDataPacked)
        {
          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Entering UAMRO");
          long long int offset = 0;

          if (hemelb::multiscale::mpwide::mpwide_comm_proc)
          {

            for (ContentsType::iterator intercommunicandData = registeredObjects.begin();
                intercommunicandData != registeredObjects.end(); intercommunicandData++)
            {
              hemelb::multiscale::Intercommunicand &sharedObject = *intercommunicandData->first;
              //std::string &label = intercommunicandData->second.second;
              IntercommunicandTypeT &resolver = *intercommunicandData->second.first;

              for (unsigned int sharedFieldIndex = 0; sharedFieldIndex < sharedObject.Values().size();
                  sharedFieldIndex++)
              {
                long long int size = GetTypeSize(resolver.Fields()[sharedFieldIndex].second,
                                                 *sharedObject.Values()[sharedFieldIndex]);
                void *buf1 = (void *) & (ICandRecvDataPacked[offset]);
                void *buf2 = (void *) & (*sharedObject.Values()[sharedFieldIndex]);

                memcpy(buf2, buf1, size);
                offset += size;
                /* Temporary diagnostic DEBUG */
                int cp = 0;
                if (hemelb::multiscale::mpwide::mpwide_comm_proc)
                {
                  cp = 1;
                }
                std::cout << "Shared value [" << sharedFieldIndex << "] = "
                    << static_cast<double*>(buf2)[sharedFieldIndex] << std::endl;
              }
            }
          }

          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("PreMPI_Bcast...%d",
                                                                                             registeredObjects.size());

          /* Over all Processes. */
          for (ContentsType::iterator intercommunicandData = registeredObjects.begin();
              intercommunicandData != registeredObjects.end(); intercommunicandData++)
          {
            hemelb::multiscale::Intercommunicand &sharedObject = *intercommunicandData->first;
            double sharedPressure = static_cast<double*>((void *) & (*sharedObject.Values()[0]))[0];
            IntercommunicandTypeT &resolver = *intercommunicandData->second.first;

            if (hemelb::multiscale::mpwide::mpwide_comm_proc)
            {
              hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Broadcasting...%lf %d",
                                                                                   sharedPressure, registeredObjects.size());
            }
            else
            {
              hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Receiving...%lf %d", sharedPressure, registeredObjects.size());
            }
            MPI_Bcast(&sharedPressure, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            static_cast<double*>((void *) & (*sharedObject.Values()[0]))[0] = sharedPressure;
            if (hemelb::multiscale::mpwide::mpwide_comm_proc)
            {
              hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Broadcasting2:...%lf",
                                                                                   sharedPressure);
            }
            else
            {
              hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Receiving2:...%lf", sharedPressure);
            }
          }

        }

        bool ExchangeWithMultiscale()
        {
          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Entering EWM");

          // 1. Pack/Serialize local shared data.
          PackRegisteredObjects(ICandSendDataPacked, registeredObjects);

          // 2. Exchange serialized shared data.
          ExchangePackages(ICandSendDataPacked, ICandRecvDataPacked);
          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Calling UAMRO");
          // 3. Unpack and merged the two serialized shared data copies.
          UnpackAndMergeRegisteredObjects(registeredObjects, ICandRecvDataPacked);

          return true;
        }

        int64_t GetRegisteredObjectsSize(ContentsType registeredObjects)
        {
          int64_t size = 0;
          int count = 0;

          for (ContentsType::iterator intercommunicandData = registeredObjects.begin();
              intercommunicandData != registeredObjects.end(); intercommunicandData++)
          {
            hemelb::multiscale::Intercommunicand &sharedObject = *intercommunicandData->first;

            //std::cout << "Number of registered objects is: " << sharedObject.Values().size() << std::endl;
            //std::string &label = intercommunicandData->second.second;
            IntercommunicandTypeT &resolver = *intercommunicandData->second.first;

            for (unsigned int sharedFieldIndex = 0; sharedFieldIndex < sharedObject.Values().size(); sharedFieldIndex++)
            {
              size += GetTypeSize(resolver.Fields()[sharedFieldIndex].second, *sharedObject.Values()[sharedFieldIndex]);

            }
            count++;
          }

          std::cout << "size obtained is: " << size << " (for " << count << " objects)." << std::endl;

          return size;
        }

        int64_t ExchangeICandDataSize(int64_t send_icand_data_size)
        {
          int64_t rsize = 0;
          if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
          {
            //std::cout << "BEFORE EXCHANGE, icand sizes are: " << send_icand_data_size << "/" << rsize << std::endl;

            MPW_SendRecv( ((char *) &send_icand_data_size),
                         sizeof(int64_t),
                         ((char *) &rsize),
                         sizeof(int64_t),
                         channels,
                         1);

            //std::cout << "EXCHANGE, icand sizes are: " << send_icand_data_size << "/" << rsize << std::endl;
          }

          return rsize;
        }

        std::map<std::string, double> &doubleContents;
        double currentTime;
        std::map<std::string, bool> & orchestration;

        /* MPWide specific parameters */
        std::string *hosts;
        int *server_side_ports;
        int *channels;
        int num_channels;
    };
  }
}

#endif // HEMELB_MULTISCALE_MPWIDE_MPWIDEINTERCOMMUNICATOR_H
