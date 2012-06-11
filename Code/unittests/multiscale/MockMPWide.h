#ifndef HEMELB_UNITTESTS_MULTISCALE_MOCKMPWIDE_H
#define HEMELB_UNITTESTS_MULTISCALE_MOCKMPWIDE_H
#include <iostream>
#include <cassert>

char* MPW_DNSResolve(char* host) {
  return "localhost";
}

/* Print all connections. */
void MPW_Print() {
  printf("MPW_Print(): MockMPWide loves connecting!");
}

/* Print the number of available streams. */
int MPW_NumChannels() 
{
  return 2;
}

/* Initialize the Cosmogrid library. */
void MPW_Init(std::string *url, int *server_side_ports, int num_channels) 
{ //this call omits client-side port binding.
  printf("MPW_Init(): last url = %s\n", url[num_channels-1].c_str()); 
  printf("MPW_Init(): last port = %d\n", server_side_ports[num_channels-1]);
}

/* Close all sockets and free data structures related to the library. */
int MPW_Finalize() {
  printf("MPW_Finalize(): MPWide finalized.\n");
}

/* Exchanges buffers between the two machines. */
void MPW_SendRecv ( char *sendbuf, long long int  sendsize, char *recvbuf, long long int  recvsize, int *channel, int num_channels) {
    if(sendsize != recvsize) {
      printf("ERROR: send and recv size mismatch in MockMPWide SendRecv: %lld, %lld!\n", sendsize, recvsize);
    }
    else {
      memcpy(recvbuf, sendbuf, sendsize); //very simple propagation of identical values.
    }
}
#endif
