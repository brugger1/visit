// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include <visit-config.h>
#include "SocketBridge.h"

#if defined(_WIN32)
#include <win32commhelpers.h>
#include <io.h>
#else
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <netdb.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <sys/socket.h>
#include <pwd.h>
#include <unistd.h>
#endif


#include <DebugStream.h>
#include <CouldNotConnectException.h>

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_BUFFER_SIZE 10000

//#define VISIT_DEBUG_SOCKETS
#ifdef VISIT_DEBUG_SOCKETS
#define DEBUG_SOCKETS(CODE) if(log != NULL){CODE ; fflush(log); }
#define DEBUG_SOCKETS_CODE(CODE) CODE
#else
#define DEBUG_SOCKETS(CODE)
#define DEBUG_SOCKETS_CODE(CODE)
#endif

static int  Listen(int port, struct sockaddr_in &sin);
static int  Accept(int listenSock, struct sockaddr_in sin);
static int  Connect(const char *host, int port);
static bool ForwardData(FILE *log, int read_fd, int write_fd, char *buff, int buffSize);
static void CloseSocket(int fd);


// ----------------------------------------------------------------------------
//                               SocketBridge
// ----------------------------------------------------------------------------

// ****************************************************************************
//  Constructor:  SocketBridge::SocketBridge
//
//  Arguments:
//    from       the new local port to listen on
//    to         the new local port to forward to
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
//  Modifications:
//    Thomas R. Treadway, Mon Oct  8 13:27:42 PDT 2007
//    Backing out SSH tunneling on Panther (MacOS X 10.3)
//
//    Gunther H. Weber, Thu Jan 14 11:38:27 PST 2010
//    Added ability to connect bridge to other host than localhost.
//
//    Brad Whitlock, Thu Oct 16 10:28:10 PDT 2014
//    Added log, logging.
//
// ****************************************************************************
SocketBridge::SocketBridge(int from, int to, const char* toHost)
{
    num_bridges = 0;
    from_port = from;
    to_port = to;
    to_host = toHost;
    buff = new char[MAX_BUFFER_SIZE];
    buffSize = MAX_BUFFER_SIZE;

    log = NULL;
    logging = false;

    listen_fd = Listen(from_port, listen_sock);
    if (listen_fd<0)
        EXCEPTION0(CouldNotConnectException);
}

// ****************************************************************************
//  Destructor:  SocketBridge::~SocketBridge
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
SocketBridge::~SocketBridge()
{
    CloseSocket(listen_fd);
    for (int i=0; i<num_bridges; i++)
    {
        CloseSocket(originating_fd[i]);
        CloseSocket(terminating_fd[i]);
    }
    delete [] buff;
    if(log != NULL)
        fclose(log);
}

// ****************************************************************************
// Method: SocketBridge::SetLogging
//
// Purpose:
//   Set whether we're doing logging or not. If we are, start a log.
//
// Arguments:
//   val : True if we're logging; false otherwise.
//
// Returns:    
//
// Note:       
//
// Programmer: Brad Whitlock
// Creation:   Thu Oct 16 10:27:15 PDT 2014
//
// Modifications:
//
// ****************************************************************************

void
SocketBridge::SetLogging(bool val)
{
    if(val)
    {
        DEBUG_SOCKETS_CODE(
        if(log == NULL)
        {
            char filename[100];
            snprintf(filename, 100, "SocketBridge_%d_%d.vlog", from_port, to_port);
            log = fopen(filename, "wt");
        }
        )
    }
    else
    {
        if(log != NULL)
        {
            fclose(NULL);
            log = NULL;
        }
    }
}

// ****************************************************************************
// Method: SocketBridge::SetBufferSize
//
// Purpose:
//   Set the buffer size.
//
// Arguments:
//   sz : The new buffer size.
//
// Returns:    
//
// Note:       This can come in handy for fixed buffer socket mode so we can
//             set the buffer size to the same as 
//             SocketConnection::FIXED_BUFFER_SIZE.
//
// Programmer: Brad Whitlock
// Creation:   Tue Oct 14 17:04:36 PDT 2014
//
// Modifications:
//
// ****************************************************************************

void
SocketBridge::SetBufferSize(int sz)
{
    if(sz > buffSize)
    {
        delete [] buff;
        buff = new char[sz];
    }
    buffSize = sz;

    DEBUG_SOCKETS(fprintf(log, "SocketBridge setting buffer size: %d.\n", sz);)
}

// ****************************************************************************
//  Method:  SocketBridge::NumActiveBridges
//
//  Purpose:
//    Return the number of active bridges.
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
int
SocketBridge::NumActiveBridges()
{
    return num_bridges;
}

// ****************************************************************************
//  Method:  SocketBridge::GetListenActivity
//
//  Purpose:
//    Returns true if the listen socket has activity.
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
bool
SocketBridge::GetListenActivity()
{
    return FD_ISSET(listen_fd, &activity);
}

// ****************************************************************************
//  Method:  SocketBridge::GetOriginatingActivity
//
//  Purpose:
//    Returns the index of the bridge that had activity on the originating
//    side (or -1 if nothing).
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
int
SocketBridge::GetOriginatingActivity()
{
    for (int i=0; i<num_bridges; i++)
    {
        if (FD_ISSET(originating_fd[i], &activity))
        {
            return i; 
        }
    }
    return -1;
}

// ****************************************************************************
//  Method:  SocketBridge::GetOriginatingActivity
//
//  Purpose:
//    Returns the index of the bridge that had activity on the terminating
//    side (or -1 if nothing).
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
int
SocketBridge::GetTerminatingActivity()
{
    for (int i=0; i<num_bridges; i++)
    {
        if (FD_ISSET(terminating_fd[i], &activity))
        {
            return i; 
        }
    }
    return -1;
}

// ****************************************************************************
//  Method:  SocketBridge::WaitForActivity
//
//  Purpose:
//    Waits for activity on either the listen socket or one of the bridges.
//    Stores the activity array for later use.  Blocking.
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
void
SocketBridge::WaitForActivity()
{
    FD_ZERO(&activity);
    int maxDescriptor = listen_fd;
    FD_SET(listen_fd, &activity);
    for (int i=0; i<num_bridges; i++)
    {
        if (originating_fd[i] > maxDescriptor)
            maxDescriptor = originating_fd[i];
        FD_SET(originating_fd[i], &activity);
    }
    for (int i=0; i<num_bridges; i++)
    {
        if (terminating_fd[i] > maxDescriptor)
            maxDescriptor = terminating_fd[i];
        FD_SET(terminating_fd[i], &activity);
    }

    int status = select(maxDescriptor+1,
                        &activity, (fd_set*)NULL, (fd_set*)NULL, NULL);
    if (status<0)
        EXCEPTION0(CouldNotConnectException);
}

// ****************************************************************************
//  Method:  SocketBridge::StartNewBridge
//
//  Purpose:
//    Start a new bridge by accepting a new incoming connection and making
//    a new outbound connection.  This assumes there is already an incoming
//    connection attempt pending.
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
//  Modifications:
//    Gunther H. Weber, Thu Jan 14 11:38:27 PST 2010
//    Added ability to connect bridge to other host than localhost.
//
// ****************************************************************************
void
SocketBridge::StartNewBridge()
{
    int ofd = Accept(listen_fd, listen_sock);
    if (ofd == -1)
        EXCEPTION0(CouldNotConnectException);
    originating_fd[num_bridges] = ofd;

    int tfd = Connect(to_host, to_port);
    if (tfd == -1)
        EXCEPTION0(CouldNotConnectException);
    terminating_fd[num_bridges] = tfd;

    DEBUG_SOCKETS(fprintf(log, "SocketBridge started new bridge.\n");)

    num_bridges++;
}

// ****************************************************************************
//  Method:  SocketBridge::CloseBridge
//
//  Purpose:
//    Close one of the bridges.  This invalidates all other bridge
//    indexes.
//
//  Arguments:
//    index      the index of the bridge to close
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
void
SocketBridge::CloseBridge(int index)
{
    CloseSocket(originating_fd[index]);
    CloseSocket(terminating_fd[index]);

    // shift the rest down
    for (int i=index; i<num_bridges-1; i++)
    {
        originating_fd[i] = originating_fd[i+1];
        terminating_fd[i] = terminating_fd[i+1];
    }

    num_bridges--;
}

// ****************************************************************************
//  Method:  SocketBridge::ForwardOrigToTerm
//
//  Purpose:
//    Forward some data from the originating end of a bridge to the
//    terminating end.
//
//  Arguments:
//    index      the index of the bridge to forward data through
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
//  Modifications:
//    Brad Whitlock, Tue Oct 14 17:00:10 PDT 2014
//    Pass in the buffer and size we'll use.
//
// ****************************************************************************
void
SocketBridge::ForwardOrigToTerm(int index)
{
    DEBUG_SOCKETS(fprintf(log, "SocketBridge::ForwardOrigToTerm\n");)
    bool success = ForwardData(log, 
                               originating_fd[index],
                               terminating_fd[index],
                               buff,
                               buffSize);
    if (!success)
        CloseBridge(index);
}

// ****************************************************************************
//  Method:  SocketBridge::ForwardTermToOrig
//
//  Purpose:
//    Forward some data from the terminating end of a bridge to the
//    originating end.
//
//  Arguments:
//    index      the index of the bridge to forward data through
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
//  Modifications:
//    Brad Whitlock, Tue Oct 14 17:00:10 PDT 2014
//    Pass in the buffer and size we'll use.
//
// ****************************************************************************
void
SocketBridge::ForwardTermToOrig(int index)
{
    DEBUG_SOCKETS(fprintf(log, "SocketBridge::ForwardTermToOrig\n");)
    bool success = ForwardData(log,
                               terminating_fd[index],
                               originating_fd[index],
                               buff,
                               buffSize);
    if (!success)
        CloseBridge(index);
}

// ****************************************************************************
//  Method:  SocketBridge::Bridge
//
//  Purpose:
//    Main loop to bridge two ports on localhost.
//
//  Arguments:
//    none
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
void
SocketBridge::Bridge()
{
    TRY
    {
        bool waiting_for_first_connection = true;
        while (NumActiveBridges() > 0 || waiting_for_first_connection)
        {
            WaitForActivity();
            if (GetListenActivity())
            {
                StartNewBridge();
                waiting_for_first_connection = false;
            }
            else
            {
                int index;
                if ((index = GetOriginatingActivity()) >= 0)
                {
                    ForwardOrigToTerm(index);
                }
                else if ((index = GetTerminatingActivity()) >= 0)
                {
                    ForwardTermToOrig(index);
                }
                else
                {
                    break;
                }
            }
        }
    }
    CATCH(CouldNotConnectException)
    {
        // Something failed, but there's not much we can do.
    }
    ENDTRY
}





// ----------------------------------------------------------------------------
//                         Socket Helper Functions
// ----------------------------------------------------------------------------




// ****************************************************************************
//  Method:  Listen
//
//  Purpose:
//    listen on a port, return the fd and fill out the sockaddr_in struct.
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
static int
Listen(int port, struct sockaddr_in &sin)
{
    int  on = 1;

    // Open a socket.
    int listenSock = socket(AF_INET, SOCK_STREAM, 0);
    if (listenSock < 0)
        return -1;

    sin.sin_family = AF_INET;
    sin.sin_addr.s_addr = htonl(INADDR_ANY);
    sin.sin_port = htons(port);
#if !defined(_WIN32)
    setsockopt(listenSock, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on));
#endif
    if (bind(listenSock, (struct sockaddr *)&sin, sizeof(sin)) < 0)
        return -1;

    if (listen(listenSock, 5) < 0)
        return -1;

    return listenSock;
}


// ****************************************************************************
//  Method:  Accept
//
//  Purpose:
//    given a listed fd and sockaddr_in struct, accept a new connection
//    and return the new fd
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
static int
Accept(int listenSock, struct sockaddr_in sin)
{
    int acceptedSock = -1;

    do
    {
#ifdef HAVE_SOCKLEN_T
        socklen_t len;
#else
        int len;
#endif
        len = sizeof(struct sockaddr);
        acceptedSock = accept(listenSock, (struct sockaddr*)&sin,&len);
    }
    while (acceptedSock == -1 && errno == EINTR);

    return acceptedSock;
}


// ****************************************************************************
//  Method:  Connect
//
//  Purpose:
//    given a host and a port, connect to it and return the new fd
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
static int
Connect(const char *host, int port)
{
#if defined(_WIN32)
    // Create a copy of the hostent struct in case other WinSock
    // functions want to modify it.
    void *hostInfo = (void *)CopyHostent(gethostbyname(host));
#else
    void *hostInfo = (void *)gethostbyname(host);
#endif

    int                connectedSock;
    struct hostent     *hp;
    struct sockaddr_in server;

    //
    // Set up the structures for opening the sockets.
    //
    hp = (struct hostent *)hostInfo;
    if (hp == NULL)
        return -1;
    memset(&server, 0, sizeof(server));
    memcpy(&(server.sin_addr), hp->h_addr, hp->h_length);
    server.sin_family = hp->h_addrtype;
    server.sin_port = htons(port);
    
    // 
    // Create a socket.
    // 
#if defined(_WIN32)
    connectedSock = socket(AF_INET, SOCK_STREAM, 0);
    if (connectedSock == INVALID_SOCKET)
    {
        return -1;
    }

    // Disable the Nagle algorithm 
    int opt = 1;
    setsockopt(connectedSock, IPPROTO_TCP, TCP_NODELAY,
               (const char FAR *)&opt, sizeof(int));

    if (connect(connectedSock, (struct sockaddr *)&server, sizeof(server)) < 0)
    {
        closesocket(connectedSock);
        return -1;
    }
#else
    connectedSock = socket(AF_INET, SOCK_STREAM, 0);
    if (connectedSock < 0)
    {
        return -1;
    }

    // Disable the Nagle algorithm 
    int opt = 1;
    setsockopt(connectedSock, IPPROTO_TCP, TCP_NODELAY, &opt, sizeof(int));
    if (connect(connectedSock, (struct sockaddr *)&server, sizeof(server)) < 0)
    {
        close(connectedSock);
        return -1;
    }
#endif

    return connectedSock;
}


// ****************************************************************************
//  Method:  ForwardData
//
//  Purpose:
//    Read a chunk of data from read_fd and send it all to write_fd.
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
//  Modification:
//    Brad Whitlock, Tue Oct 14 16:56:21 PDT 2014
//    Pass in the buffer to use and its size. Added logging.
//
// ****************************************************************************
static bool
ForwardData(FILE *log, int read_fd, int write_fd, char *buff, int buffSize)
{
    DEBUG_SOCKETS_CODE(const char *mName = "ForwardData: ";)
    DEBUG_SOCKETS(fprintf(log, "%s begin\n", mName);
            fprintf(log, "%sread(%d): buffSize=%d\n", mName, read_fd, buffSize);)
#ifndef WIN32
    int nread = read(read_fd, buff, buffSize);
#else
    int nread = _read(read_fd, buff, buffSize);
#endif
    DEBUG_SOCKETS(fprintf(log, "%sread(%d) returned %d bytes ", mName, read_fd, nread);)
    
    if (nread <= 0)
        return false;

    DEBUG_SOCKETS(
        int n = std::min(nread, 10);
        fprintf(log, "{");
        for(int i = 0; i < n; ++i)
            fprintf(log, "%d, ", int(buff[i]));
        if(n < nread)
            fprintf(log, "...");
        fprintf(log, "}\n");
    )

    size_t          nleft, nwritten;
    const char      *ptr;
    ptr = (const char*)buff;
    nleft = nread;
    while (nleft > 0)
    {
        DEBUG_SOCKETS(fprintf(log, "%ssend(%d) buffSize=%d\n", mName, write_fd, buffSize);)

        if ((nwritten = send(write_fd, ptr, nleft, 0)) <= 0)
            break;

        DEBUG_SOCKETS(fprintf(log, "%ssend(%d) sent %d bytes\n", mName, write_fd, (int)nwritten);)

        nleft -= nwritten;
        ptr   += nwritten;
    }

    if (nleft > 0)
    {
        DEBUG_SOCKETS(fprintf(log, "%serror\n", mName);)
        return false;
    }

    DEBUG_SOCKETS(fprintf(log, "%send\n", mName);)

    return true;
}


// ****************************************************************************
//  Method:  CloseSocket
//
//  Purpose:
//    Close a socket file descriptor.
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 24, 2007
//
// ****************************************************************************
static void
CloseSocket(int fd)
{
    if (fd < 0)
        return;

#if defined(_WIN32)
    closesocket(fd);
#else
    close(fd);
#endif
}

