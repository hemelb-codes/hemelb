#ifndef HEME_STEERING_HTTPPOST_H
#define HEME_STEERING_HTTPPOST_H

namespace heme
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

#endif // HEME_STEERING_HTTPPOST_H
