#ifndef HTTP_POST
#define HTTP_POST

class HTTP
{
  public:
    // Send an HTTP request
    static int request (char* hostname, short port, char* api, char* resourceid, char* parameters);

    // Get host details
    static void get_host_details(char* rank_0_host_details, char* ip_addr);
};

#endif
