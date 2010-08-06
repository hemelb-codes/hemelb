#ifndef __steering_HttpPost_h_
#define __steering_HttpPost_h_

class HttpPost
{
 public:
  // Send an HTTP request
  static int request (char* hostname, short port,
		      char* api, char* resourceid, char* parameters);
  
  // Get host details
  static void get_host_details(char* rank_0_host_details, char* ip_addr);
};

#endif//__steering_HttpPost_h_
