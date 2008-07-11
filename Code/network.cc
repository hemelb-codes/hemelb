#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>

#include <sys/stat.h>

int recv_all (int sockid, char *buf, int *length) {

	int received_bytes = 0;
	int bytes_left_to_receive = *length;
	int n;

	while (received_bytes < *length) {
		n = recv(sockid, buf+received_bytes, bytes_left_to_receive, 0);
		if (n == -1) break;
		received_bytes += n;
		bytes_left_to_receive -= n;
	}

	*length = received_bytes;
	return n == -1 ? -1 : 0;
}


int send_all(int sockid, char *buf, int *length ) {
  
	int sent_bytes = 0;
	int bytes_left_to_send = *length;
	int n;
	
	while( sent_bytes < *length ) {
		n = send(sockid, buf+sent_bytes, bytes_left_to_send, 0);
		if (n == -1) break;
		sent_bytes += n;
		bytes_left_to_send -= n;
	}
	
	*length = sent_bytes;
	
	return n== -1 ? -1 : 0;
}

