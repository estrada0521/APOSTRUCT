/* TCP sockets for the GUI server and the show client, POSIX or Winsock. */

#ifndef APOSTRUCT_NET_H
#define APOSTRUCT_NET_H

#ifdef _WIN32
#include <winsock2.h>
#include <ws2tcpip.h>

typedef SOCKET net_socket;
#define NET_INVALID INVALID_SOCKET
static inline int net_init(void)
{
    WSADATA data;
    return WSAStartup(MAKEWORD(2, 2), &data);
}
static inline int net_close(net_socket s) { return closesocket(s); }
#else
#include <arpa/inet.h>
#include <netdb.h>
#include <netinet/in.h>
#include <signal.h>
#include <sys/socket.h>
#include <unistd.h>

typedef int net_socket;
#define NET_INVALID (-1)
/* a peer that hangs up must not end the process */
static inline int net_init(void)
{
    signal(SIGPIPE, SIG_IGN);
    return 0;
}
static inline int net_close(net_socket s) { return close(s); }
#endif

#endif
