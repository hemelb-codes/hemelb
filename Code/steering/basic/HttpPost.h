// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_STEERING_BASIC_HTTPPOST_H
#define HEMELB_STEERING_BASIC_HTTPPOST_H

namespace hemelb
{
  namespace steering
  {

    class HttpPost
    {
      public:
        /**
         * Send an HTTP Request.
         *
         * @param hostname The name of the host.
         * @param port The port to send to.
         * @param api The API name to include in the request.
         * @param resourceid An identifier for the desired resource.
         * @return Error code. 0 on success.
         */
        static int request(const char* hostname,
                           const in_port_t port,
                           const char* api,
                           const char* resourceid);

      private:
        /**
         *Get details and ip address for this host.
         *
         * @param rank_0_host_details
         * @param ip_addr
         */
        static void get_host_details(char* rank_0_host_details);

        /**
         * Helper method to send a request.
         *
         * @param iSocket
         * @param iMessage
         * @return
         */
        static ssize_t Send_Request(int iSocket, const char *iMessage);

        static const int MYPORT = 65250;
    };

  }
}

#endif // HEMELB_STEERING_BASIC_HTTPPOST_H
