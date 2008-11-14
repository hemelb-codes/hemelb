#define _DEBUG_PRINT(X)   /* X */

#include <iostream>
#include <string>
#include <cstdlib>
#include <cassert>
#include <netdb.h>
#include <cstring>


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

#define MYPORT 65250

using namespace std;


void get_host_details(char* rank_0_host_details, char* ip_addr) {

#ifdef NGS2Leeds
	char hostname[256];
	gethostname(hostname, 256);
        sprintf(rank_0_host_details, "%s.ngs.leeds.ac.uk:%i", hostname, MYPORT);
	fprintf(stderr, "MPI rank 0 public interface details - %s\n", rank_0_host_details);
#elif NGS2Manchester
	char hostname[256];
	gethostname(hostname, 256);
        sprintf(rank_0_host_details, "%s.vidar.ngs.manchester.ac.uk:%i", hostname, MYPORT);
	fprintf(stderr, "MPI rank 0 public interface details - %s\n", rank_0_host_details);
#elif CCS
	char hostname[256];
	gethostname(hostname, 256);
        sprintf(rank_0_host_details, "%s.chem.ucl.ac.uk:%i", hostname, MYPORT);
	fprintf(stderr, "MPI rank 0 public interface details - %s\n", rank_0_host_details);
#elif LONI
	char hostname[256];
	gethostname(hostname, 256);
        sprintf(rank_0_host_details, "%s:%i", hostname, MYPORT);
	fprintf(stderr, "MPI rank 0 public interface details - %s\n", rank_0_host_details);
#else

        struct utsname  name;
        struct hostent *host;

        if(uname(&name) < 0) {
                fprintf(stderr, "STEER: Get_fully_qualified_hostname: uname failed\n");
                exit(0x0);
        }

        host = gethostbyname(name.nodename);

	for(int address_id = 0; address_id < 4  ; address_id++ ) {

	if( host->h_addr_list[address_id] == NULL ) break;

		printf("checking Address ID %i..\n", address_id);

        if(host->h_length == 4){

                sprintf(ip_addr, "%d.%d.%d.%d",
	                (unsigned char)(host->h_addr_list[address_id][0]),
   	             (unsigned char)(host->h_addr_list[address_id][1]),
   	             (unsigned char)(host->h_addr_list[address_id][2]),
   	             (unsigned char)(host->h_addr_list[address_id][3]));

				printf("NAME %s IP %s\n", host->h_name, ip_addr);

				if( (unsigned char)(host->h_addr_list[address_id][0]) == 10 ) {
					printf("IP %s is not public..\n", ip_addr);
					continue;
				}

				if( (unsigned char)(host->h_addr_list[address_id][0]) == 127 ) {
					printf("IP %s is not public..\n", ip_addr);
					continue;
				}

				if( ( (unsigned char)(host->h_addr_list[address_id][0]) == 172 ) &&
				    ( (unsigned char)(host->h_addr_list[address_id][1]) == 16 ) ) {
					printf("IP %s is not public..\n", ip_addr);
					continue;
				}

				if( ( (unsigned char)(host->h_addr_list[address_id][0]) == 192 ) &&
				    ( (unsigned char)(host->h_addr_list[address_id][1]) == 168 ) &&
				    ( (unsigned char)(host->h_addr_list[address_id][1]) == 0  ) ) {
					printf("IP %s is not public..\n", ip_addr);
					continue;
				}

        } else {

				printf("address is not IP4..\n");
				continue;

	}
	}

        printf("Public hostname to use will be %s IP %s\n", host->h_name, ip_addr);

        sprintf(rank_0_host_details, "%s:%i", host->h_name, MYPORT);

#endif

}

/* It's OK, this was written by a Java developer */

#define SEND_RQ(MSG) \
                /*cout<<send_str;*/ \
  send(sock,MSG,strlen(MSG),0);

int request (char* hostname, short port, char* api, char* resourceid, char* parameters) {

	sockaddr_in sin;
	int sock = socket (AF_INET, SOCK_STREAM, 0);
	if (sock == -1) {
		return -100;
	}

    sin.sin_family = AF_INET;
    sin.sin_port = htons((unsigned short)port);

    struct hostent * host_addr = gethostbyname(hostname);

    if(host_addr==NULL) {
      _DEBUG_PRINT( cout<<"Unable to locate host"<<endl );
      return -103;
    }

    sin.sin_addr.s_addr = *((int*)*host_addr->h_addr_list) ;
    _DEBUG_PRINT( cout<<"Port :"<<sin.sin_port<<", Address : "<< sin.sin_addr.s_addr<<endl);

    if( connect (sock,(const struct sockaddr *)&sin, sizeof(sockaddr_in) ) == -1 ) {
     _DEBUG_PRINT( cout<<"connect failed"<<endl ) ;
     return -101;
    }

 string send_str;

 SEND_RQ("POST ");
 SEND_RQ(api);
 SEND_RQ(resourceid);
 SEND_RQ(" HTTP/1.0\r\n");
 SEND_RQ("Accept: */*\r\n");
 SEND_RQ("User-Agent: Mozilla/4.0\r\n");

 char content_header[100];
 sprintf(content_header,"Content-Length: %d\r\n",strlen(parameters));
 SEND_RQ(content_header);
 SEND_RQ("Accept-Language: en-us\r\n");
 SEND_RQ("Accept-Encoding: gzip, deflate\r\n");
 SEND_RQ("Host: ");
 SEND_RQ("hostname");
 SEND_RQ("\r\n");
 SEND_RQ("Content-Type: text/plain\r\n");
 
 //If you need to send a basic authorization
 //string Auth        = "username:password";
 //Figureout a way to encode test into base64 !
 //string AuthInfo    = base64_encode(reinterpret_cast<const unsigned char*>(Auth.c_str()),Auth.length());
 //string sPassReq    = "Authorization: Basic " + AuthInfo;
 //SEND_RQ(sPassReq.c_str());

 SEND_RQ("\r\n");
 // SEND_RQ("\r\n");
 SEND_RQ(parameters);
 SEND_RQ("\r\n");

 _DEBUG_PRINT(cout<<"####HEADER####"<<endl);
 char c1[1];
 int l,line_length;
 bool loop = true;
 bool bHeader = false;

 string message;

 while(loop) {
   l = recv(sock, c1, 1, 0);
   if(l<0) loop = false;
   if(c1[0]=='\n') {
       if(line_length == 0) loop = false;

       line_length = 0;
       if(message.find("201 Created") != string::npos)
	       bHeader = true;

   }
   else if(c1[0]!='\r') line_length++;
   _DEBUG_PRINT( cout<<c1[0]);
   message += c1[0];
 }

 message="";
 if(bHeader) {

     _DEBUG_PRINT( cout<<"####BODY####"<<endl) ;
     char p[1024];
     while((l = recv(sock,p,1023,0)) > 0)  {
         _DEBUG_PRINT( cout.write(p,l)) ;
	     p[l] = '\0';
	     message += p;
     }

     _DEBUG_PRINT( cout << message.c_str());
 } else {
	 return -102;
 }

 return 0;

}

