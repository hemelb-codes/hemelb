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
      static int request (char* hostname, short port,
			  char* api, char* resourceid, char* parameters);
  
      // Get host details
      static void get_host_details(char* rank_0_host_details, char* ip_addr);
    };

  }
}

#endif // HEMELB_STEERING_ON_HTTPPOST_H
