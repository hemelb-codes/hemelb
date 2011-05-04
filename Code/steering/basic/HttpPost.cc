#include <iostream>
#include <string>
#include <cstdlib>
#include <cassert>
#include <netdb.h>
#include <cstring>
#include <cstdio>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>
#include <netdb.h>
#include <sys/utsname.h>

#include "log/Logger.h"
#include "steering/basic/HttpPost.h"

using namespace std;

namespace hemelb
{
  namespace steering
  {

    void HttpPost::get_host_details(char* rank_0_host_details)
    {
      // This really does modify the strings passed in. Without
      // checking their lengths, of course.
      char hostname[256];
      char ip_addr[16];

      // On specific machines, just get the host name and insert it into the rank_0_host_details parameter
#ifdef NGS2Leeds
      gethostname(hostname, 256);
      sprintf(rank_0_host_details, "%s.ngs.leeds.ac.uk:%i", hostname, MYPORT);
#elif NGS2Manchester
      gethostname(hostname, 256);
      sprintf(rank_0_host_details, "%s.vidar.ngs.manchester.ac.uk:%i", hostname, MYPORT);
#elif CCS
      gethostname(hostname, 256);
      std::sprintf(rank_0_host_details, "%s.chem.ucl.ac.uk:%i", hostname, MYPORT);
#elif LONI
      gethostname(hostname, 256);
      sprintf(rank_0_host_details, "%s:%i", hostname, MYPORT);
#elif NCSA
      gethostname(hostname, 256);
      sprintf(rank_0_host_details, "%s.ncsa.teragrid.org:%i", hostname, MYPORT);
#else

      // If not on a specific known machine, need to do something more clever.
      struct utsname name;

      if (uname(&name) < 0)
      {
        fprintf(stderr, "STEER: Get_fully_qualified_hostname: uname failed\n");
        exit(0x0);
      }

      // Get information about our host.
      struct hostent *host = gethostbyname(name.nodename);

      // Now go through every associated IP adress, and use the first valid, public one.
      for (int address_id = 0; address_id < 4; address_id++)
      {
        // If no address, break.
        if (host->h_addr_list[address_id] == NULL)
          break;

        printf("checking Address ID %i...\n", address_id);

        // If IP4, print details.
        if (host->h_length != 4)
        {
          printf("address is not IP4..\n");
          continue;
        }

        sprintf(ip_addr,
                "%d.%d.%d.%d",
                (unsigned char) (host->h_addr_list[address_id][0]),
                (unsigned char) (host->h_addr_list[address_id][1]),
                (unsigned char) (host->h_addr_list[address_id][2]),
                (unsigned char) (host->h_addr_list[address_id][3]));

        printf("NAME %s IP %s\n", host->h_name, ip_addr);

        // Private addresses (see RFC-1918).
        if ((unsigned char) (host->h_addr_list[address_id][0]) == 10)
        {
          printf("IP %s is not public..\n", ip_addr);
          continue;
        }

        // Loopback addresses.
        if ((unsigned char) (host->h_addr_list[address_id][0]) == 127)
        {
          printf("IP %s is not public..\n", ip_addr);
          continue;
        }

        // Private addresses.
        if ((unsigned char) (host->h_addr_list[address_id][0]) == 172
            && (unsigned char) (host->h_addr_list[address_id][1]) >= 16
            && (unsigned char) (host->h_addr_list[address_id][1]) < 32)
        {
          printf("IP %s is not public..\n", ip_addr);
          continue;
        }

        // Private IP addresses.
        if ((unsigned char) (host->h_addr_list[address_id][0]) == 192
            && (unsigned char) (host->h_addr_list[address_id][1]) == 168)
        {
          printf("IP %s is not public..\n", ip_addr);
          continue;
        }
      }

      sprintf(rank_0_host_details, "%s:%i (IP %s)", host->h_name, MYPORT, ip_addr);
#endif

      log::Logger::Log<log::Info, log::Singleton>("MPI public interface details - %s",
                                                  rank_0_host_details);
    }

    ssize_t HttpPost::Send_Request(int iSocket, const char *iMessage)
    {
      return send(iSocket, iMessage, strlen(iMessage), 0);
    }

    int HttpPost::request(const char* hostname,
                          const in_port_t port,
                          const char* api,
                          const char* resourceid)
    {
      // Get the host name to communicate with.
      char host_name[1024];
      get_host_details(host_name);

      // Attempt to get a socket to use.
      int sock = socket(AF_INET, SOCK_STREAM, 0);
      if (sock == -1)
      {
        return -100;
      }

      sockaddr_in sin;
      sin.sin_family = AF_INET;
      sin.sin_port = (in_port_t) htonl(port);

      // Get name for the other end of the connection.
      struct hostent * host_addr = gethostbyname(hostname);

      if (host_addr == NULL)
      {
        return -103;
      }

      sin.sin_addr.s_addr = * ((int*) *host_addr->h_addr_list);

      // Attempt to connect.
      if (connect(sock, (const struct sockaddr *) &sin, sizeof(sockaddr_in)) == -1)
      {
        return -101;
      }

      // Now we perform the actual sending.
      Send_Request(sock, "POST ");
      Send_Request(sock, api);
      Send_Request(sock, resourceid);
      Send_Request(sock, " HTTP/1.0\r\n");
      Send_Request(sock, "Accept: */*\r\n");
      Send_Request(sock, "User-Agent: Mozilla/4.0\r\n");

      char content_header[100];
      sprintf(content_header, "Content-Length: %d\r\n", int (strlen(host_name)));

      Send_Request(sock, content_header);
      Send_Request(sock, "Accept-Language: en-us\r\n");
      Send_Request(sock, "Accept-Encoding: gzip, deflate\r\n");
      Send_Request(sock, "Host: ");
      Send_Request(sock, "hostname");
      Send_Request(sock, "\r\n");
      Send_Request(sock, "Content-Type: text/plain\r\n");
      Send_Request(sock, "\r\n");
      Send_Request(sock, host_name);
      Send_Request(sock, "\r\n");

      /* If you need to send a basic authorization
       *  string Auth        = "username:password";
       *  Figureout a way to encode test into base64 !
       *  string AuthInfo    = base64_encode(reinterpret_cast<const unsigned char*>(Auth.c_str()),Auth.length());
       *  string sPassReq    = "Authorization: Basic " + AuthInfo;
       *  Send_Request(sock, sPassReq.c_str());
       *  */

      // Receive 1 character at a time until a whole response is recovered.
      int line_length = 0;
      bool loop = true;
      bool bHeader = false;

      string message;

      while (loop)
      {
        char c1[1];

        // receive 1 char from the socket
        ssize_t l = recv(sock, c1, 1, 0);

        // An error occurred.
        if (l < 0)
          loop = false;

        if (c1[0] == '\n')
        {
          // Received a newline. If the only character on this line, we're at the end of the
          // response.
          if (line_length == 0)
            loop = false;

          line_length = 0;

          // if the message contains a success response code, we are in the header.
          if (message.find("201 Created") != string::npos)
            bHeader = true;
        }
        else if (c1[0] != '\r')
        {
          // Didn't get newline and didn't get carriage return.
          line_length++;
        }

        message += c1[0];
      }

      message = "";

      // If we are indeed reading a header, read the rest of the message.
      if (bHeader)
      {
        char p[1024];

        ssize_t l;

        while ( (l = recv(sock, p, 1023, 0)) > 0)
        {
          p[l] = '\0';
          message += p;
        }
      }
      else
      {
        return -102;
      }

      return 0;
    }

  }
}

