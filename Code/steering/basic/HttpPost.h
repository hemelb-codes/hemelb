#ifndef HEMELB_STEERING_ON_HTTPPOST_H
#define HEMELB_STEERING_ON_HTTPPOST_H

namespace hemelb
{
  namespace steering
  {

    class HttpPost
    {
      public:
        // Send an HTTP request
        static int request(const char* hostname,
                           const short port,
                           const char* api,
                           const char* resourceid,
                           const char* parameters);

        // Get host details
        static void get_host_details(char* rank_0_host_details, char* ip_addr);

      private:
        static int Send_Request(int iSocket, const char *iMessage);

        static const int MYPORT = 65250;
    };

  }
}

#endif // HEMELB_STEERING_ON_HTTPPOST_H
