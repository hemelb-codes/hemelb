#ifndef __steer_http_post_h_
#define __steer_http_post_h_

class HTTP
{
  public:
    // Send an HTTP request
    static int request (char* hostname, short port, char* api, char* resourceid, char* parameters);

    // Get host details
    static void get_host_details(char* rank_0_host_details, char* ip_addr);
};

#endif//__steer_http_post_h_
