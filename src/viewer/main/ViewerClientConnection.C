// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include <ViewerClientConnection.h>

#include <QSocketNotifier>
#include <QTimer>

#include <AttributeSubject.h>
#include <ExistingRemoteProcess.h>
#include <LostConnectionException.h>
#include <ParentProcess.h>
#include <RemoteProcess.h>
#include <Xfer.h>
#include <ViewerState.h>
#include <ViewerRPC.h>
#include <WebSocketConnection.h>

#include <DebugStream.h>

// ****************************************************************************
// Method: ViewerClientConnection::ViewerClientConnection
//
// Purpose: 
//   ViewerClientConnection constructor.
//
// Arguments:
//   s      : The viewer state object to use as a template for this connection's
//            viewer state.
//   parent : The object's parent.
//   name   : The name of the object.
//
// Programmer: Brad Whitlock
// Creation:   Tue May 31 13:36:05 PST 2005
//
// Modifications:
//   Brad Whitlock, Mon Feb 12 17:57:18 PST 2007
//   Changed base class.
//
//   Brad Whitlock, Fri May  9 14:42:50 PDT 2008
//   Added name.
//
// ****************************************************************************

ViewerClientConnection::ViewerClientConnection(const ViewerState *s,
    QObject *parent, const QString &n,const bool _allState) : ViewerBaseUI(parent),
    SimpleObserver(), name(n)
{
    notifier = 0;
    ownsNotifier = false;
    remoteProcess = 0;
    parentProcess = 0;
    initialStateStage = 0;
    allState = _allState;
    viewerState = new ViewerState(*s);

    // Hook Xfer up to the objects in the viewerState object.
    xfer = new Xfer;
    for(int i = 0; i < viewerState->GetNumStateObjects(); ++i)
    {
        xfer->Add(viewerState->GetStateObject(i));
        viewerState->GetStateObject(i)->Attach(this);
    }
    xfer->CreateNewSpecialOpcode(); // animationStopOpcode
    xfer->CreateNewSpecialOpcode(); // iconifyOpcode
    /// the local attributes are clones, so copy the guido
    clientAtts.SetGuido(viewerState->GetViewerClientAttributes()->GetGuido());
}

// ****************************************************************************
// Method: ViewerClientConnection::ViewerClientConnection
//
// Purpose: 
//   ViewerClientConnection constructor.
//
// Arguments:
//   p      : The ParentProcess object to be used by this object.
//   sn     : The socket notifier to be used by this object.
//   s      : The viewer state to be used by this object.
//   parent : The object's parent.
//   n      : The name of the object.
//
// Note:       We use this constructor to repackage some existing objects
//             that we use to set up the connection to the original client
//             as a client connection once the viewer is all set up.
//
// Programmer: Brad Whitlock
// Creation:   Tue May 31 13:38:17 PST 2005
//
// Modifications:
//   Brad Whitlock, Thu Aug 14 13:19:41 PDT 2008
//   It's possible that the notifier and the parent process will be 0 if 
//   we're using the viewer without a client.
//
// ****************************************************************************

ViewerClientConnection::ViewerClientConnection(ParentProcess *p,
    QSocketNotifier *sn, const ViewerState *s, QObject *parent,
    const QString &n, const bool _allState) : ViewerBaseUI(parent), name(n)
{
    notifier = sn;
    if(notifier != 0)
    {
        connect(notifier, SIGNAL(activated(int)),
                this, SLOT(ReadFromClientAndProcess(int)));
    }
    ownsNotifier = false;
    allState = _allState;
    remoteProcess = 0;
    parentProcess = p;
    initialStateStage = 0;

    // Hook Xfer up to the objects in the viewerState object.
    xfer = new Xfer;
    if(parentProcess != 0)
    {
        xfer->SetInputConnection(parentProcess->GetWriteConnection());
        xfer->SetOutputConnection(parentProcess->GetReadConnection());
    }
    viewerState = new ViewerState(*s);
    for(int i = 0; i < viewerState->GetNumStateObjects(); ++i)
    {    
        xfer->Add(viewerState->GetStateObject(i));
        viewerState->GetStateObject(i)->Attach(this);
    }
    xfer->CreateNewSpecialOpcode(); // animationStopOpcode
    xfer->CreateNewSpecialOpcode(); // iconifyOpcode
    /// the local attributes are clones, so copy the guido
    clientAtts.SetGuido(viewerState->GetViewerClientAttributes()->GetGuido());
}

// ****************************************************************************
// Method: ViewerClientConnection::~ViewerClientConnection
//
// Purpose: 
//   Destructor for the ViewerClientConnection class.
//
// Programmer: Brad Whitlock
// Creation:   Tue May 31 13:40:48 PST 2005
//
// Modifications:
//   Brad Whitlock, Mon Jul 11 08:43:49 PDT 2005
//   Fixed a memory problem on Windows.
//
//   Brad Whitlock, Thu Aug 14 13:19:22 PDT 2008
//   Don't assume that the notifier will exist always.
//
// ****************************************************************************

ViewerClientConnection::~ViewerClientConnection()
{
    if(ownsNotifier)
        delete notifier;

    delete xfer;
    delete viewerState;
    delete remoteProcess;
    delete parentProcess;
}

// ****************************************************************************
// Method: ViewerClientConnection::LaunchClient
//
// Purpose: 
//   This method launches the specified client and initializes the
//   ViewerClientConnection object so that it is connected to the launched
//   client.
//
// Arguments:
//   program : The name of the program to launch.
//   args    : A list of arguments to tbe passed to the launched client.
//   cb      : Connect callback function. If it is 0 then we launch a new
//             process otherwise we use an ExistingRemoteProcess object.
//   cbData  : callback data for cb.
//   connectProgressCB : The event processing callback to use while we're
//                       launching the new client.
//   connectProgressCBData : Data for the connectProgressCB callback.
//
// Note:       
//
// Programmer: Brad Whitlock
// Creation:   Tue May 31 13:41:13 PST 2005
//
// Modifications:   
//   Jeremy Meredith, Thu May 24 10:35:14 EDT 2007
//   Added SSH tunneling option to RemoteProcess::Open, and set it to false.
//
//   Brad Whitlock, Fri May 23 11:08:53 PDT 2008
//   Qt 4.
//
//   Jeremy Meredith, Thu Feb 18 15:25:27 EST 2010
//   Split HostProfile int MachineProfile and LaunchProfile.
//
//   Brad Whitlock, Thu Feb 24 23:52:23 PST 2011
//   Send state objects out 1 by 1 on Windows to make it less likely that we
//   clog the socket. It also gives more time for the gui to make requests
//   from the viewer.
//
//   Eric Brugger, Mon May  2 17:05:03 PDT 2011
//   I added the ability to use a gateway machine when connecting to a
//   remote host.
//
//   Brad Whitlock, Tue Jun  5 16:35:01 PDT 2012
//   Pass MachineProfile::Default output to RemoteProcess.
//
//   Kathleen Biagas, Thu Aug 30 2018
//   Reworked logic for sending State objects on Windows.  The 1 by 1 method
//   when launching cli from command window is a problem, so limit that path
//   to launch of gui.
//
// ****************************************************************************

void
ViewerClientConnection::LaunchClient(const std::string &program,
    const stringVector &args, 
    void (*cb)(const std::string &, const stringVector &, void *),
    void *cbData,
    bool (*connectProgressCB)(void *, int),
    void *connectProgressCBData)
{
    const char *mName = "ViewerClientConnection::LaunchClient: ";

    if(parentProcess != 0)
    {
        debug1 << mName << "This method is only to be used when "
            "listening for new client connections." << endl;
        return;
    }

    if(cb == 0)
    {
        debug1 << mName << "Creating RemoteProcess for "
               << program.c_str() << endl;
        remoteProcess = new RemoteProcess(program);
    }
    else
    {
        debug1 << mName << "Creating ExistingRemoteProcess for "
               << program.c_str() << endl;
        ExistingRemoteProcess *erp = new ExistingRemoteProcess(program);
        erp->SetConnectCallback(cb);
        erp->SetConnectCallbackData(cbData);
        remoteProcess = erp;
    }

    debug1 << mName << "Process arguments:" << endl;
    debug1 << "\t-reverse_launch" << endl;
    remoteProcess->AddArgument("-reverse_launch");
    for(size_t i = 0; i < args.size(); ++i)
    {
        remoteProcess->AddArgument(args[i]);
        debug1 << "\t" << args[i].c_str() << endl;
    }

    // Set a progress callback for the remote process object.
    if(connectProgressCB != 0)
    {
        remoteProcess->SetProgressCallback(connectProgressCB, 
                                           connectProgressCBData);
    }

    // Try opening the client.
    debug1 << mName << "Opening client connection." << endl;
    remoteProcess->Open(MachineProfile::Default(), 1, 1);

    debug1 << mName << "Successfully opened client connection." << endl;

    // Hook up the remote process's input/output to xfer.
    xfer->SetInputConnection(remoteProcess->GetWriteConnection());
    xfer->SetOutputConnection(remoteProcess->GetReadConnection());

    //
    // Create a QSocketNotifier that tells us to call ReadFromClientAndProcess.
    //
    if(remoteProcess->GetWriteConnection())
    {
        int desc = remoteProcess->GetWriteConnection()->GetDescriptor();
        if(desc != -1)
        {
            debug1 << mName << "Creating socket notifier to listen to client."
                   << endl;
            WebSocketConnection* conn = dynamic_cast<WebSocketConnection*>(remoteProcess->GetWriteConnection());

            if(conn) {
                connect(conn, SIGNAL(frameRead(int)), this, SLOT(ReadFromClientAndProcess(int)));
                connect(conn, SIGNAL(disconnected()), this, SLOT(ForceDisconnectClient()));
            }
            else {
                notifier = new QSocketNotifier(desc, QSocketNotifier::Read, 0);
                connect(notifier, SIGNAL(activated(int)),
                        this, SLOT(ReadFromClientAndProcess(int)));
                ownsNotifier = true;
            }

            /// if there is a close from client
//            QTcpSocket* s = new QTcpSocket();
//            s->setSocketDescriptor(desc);
//            connect(s, SIGNAL(aboutToClose()), SLOT(ForceDisconnectClient()));
        }
    }
    bool sendOnlyInitialState = false;
#ifdef _WIN32
    // Sending only initial state is a problem if launching cli via Command
    // windows, so limit it to only gui.  Not sure if this is truly still
    // necessary for the gui, but will leave it in nevertheless.
    if (std::find(args.begin(), args.end(), "-gui") != args.end())
        sendOnlyInitialState = true;
#endif
    if (sendOnlyInitialState)
    {
        // Initiate sending state objects to the client.
        initialStateStage = allState ? 0 : viewerState->FreelyExchangedState();
        QTimer::singleShot(50, this, SLOT(sendInitialState()));
    }
    else
    {
        // Send all of the state except for the first 7 state objects, which
        // are: ViewerRPC, PostponedRPC, syncAtts, messageAtts, statusAtts,
        // metaData, silAtts.
        debug1 << mName << "Sending state objects to client." << endl;
        for(int i = allState ? 0 : viewerState->FreelyExchangedState();
            i < viewerState->GetNumStateObjects(); ++i)
        {
            viewerState->GetStateObject(i)->SelectAll();
            SetUpdate(false);
            if(allState) viewerState->GetStateObject(i)->SetSendMetaInformation(true);
            viewerState->GetStateObject(i)->Notify();
            if(allState) viewerState->GetStateObject(i)->SetSendMetaInformation(false);
        }
    }

    debug1 << mName << "Done" << endl;
}

// ****************************************************************************
// Method: ViewerClientConnection::sendInitialState
//
// Purpose: 
//   This is a Qt slot function that we can use to send state out to the client
//   one state object at a time from the main event loop.
//
// Programmer: Brad Whitlock
// Creation:   Thu Feb 24 23:52:23 PST 2011
//
// Modifications:
//   
// ****************************************************************************

void
ViewerClientConnection::sendInitialState()
{
    // Send one state object.
    viewerState->GetStateObject(initialStateStage)->SelectAll();
    SetUpdate(false);
    if(allState) viewerState->GetStateObject(initialStateStage)->SetSendMetaInformation(true);
    viewerState->GetStateObject(initialStateStage)->Notify();
    if(allState) viewerState->GetStateObject(initialStateStage)->SetSendMetaInformation(false);

    // See if we should send another state object in a deferred manner.
    initialStateStage++;
    if(initialStateStage < viewerState->GetNumStateObjects())
        QTimer::singleShot(50, this, SLOT(sendInitialState()));
}

// ****************************************************************************
// Method: ViewerClientConnection::SetupSpecialOpcodeHandler
//
// Purpose: 
//   Sets up a special opcode handler with xfer so interrupt, stop, hide
//   all work.
//
// Arguments:
//   cb   : The special opcode handler.
//   data : The callback data for the opcode handler.
//
// Programmer: Brad Whitlock
// Creation:   Tue May 31 13:45:30 PST 2005
//
// Modifications:
//   
// ****************************************************************************

void
ViewerClientConnection::SetupSpecialOpcodeHandler(void (*cb)(int,void*),
    void *data)
{
    xfer->SetupSpecialOpcodeHandler(cb, data);
}

// ****************************************************************************
// Method: ViewerClientConnection::Update
//
// Purpose: 
//   This method gets called when a state object is read by xfer. We use the
//   opportunity to emit signals that tell the ViewerSubject to add the input
//   to its main input connection or to disconnect this client.
//
// Arguments:
//   subj : The subject that was notified.
//
// Programmer: Brad Whitlock
// Creation:   Tue May 31 13:48:44 PST 2005
//
// Modifications:
//   Brad Whitlock, Thu Aug 14 13:17:42 PDT 2008
//   Don't assume the notifier will exist.
//
// ****************************************************************************

void
ViewerClientConnection::Update(Subject *subj)
{
    bool doEmit = true;

    //
    // Before we tell the viewersubject to process this input, check to
    // see if it is a viewer RPC and if it is DetachRPC then emit a
    // different signal.
    //
    if(subj == viewerState->GetStateObject(0))
    {
        ViewerRPC *rpc = (ViewerRPC *)subj;
        if(rpc->GetRPCType() == ViewerRPC::DetachRPC)
        {
            doEmit = false;
            if(notifier != 0)
            {
                disconnect(notifier, SIGNAL(activated(int)),
                           this, SLOT(ReadFromClientAndProcess(int)));
                notifier->setEnabled(false);
            }
            emit DisconnectClient(this);
        }

        if(rpc->GetRPCType() == ViewerRPC::ExportRPC) {
            /// embed the response to the specific client..
            rpc->SetIntArg1(clientAtts.GetId());
        }
    }
    
    if(doEmit)
    {
        emit InputFromClient(this, (AttributeSubject *)subj);
    }
}

// ****************************************************************************
// Method: ViewerClientConnection::BroadcastToClient
//
// Purpose: 
//   This method broadcasts the specified state object to the client.
//
// Arguments:
//   src : The state object to send to the client.
//
// Programmer: Brad Whitlock
// Creation:   Tue May 31 13:50:14 PST 2005
//
// Modifications:
//   
// ****************************************************************************

void
ViewerClientConnection::BroadcastToClient(AttributeSubject *src)
{
    // ViewerState and xfer, which sets an objects "guido", have the
    // same indexing. The AttributeSubject passed is registered with
    // the viewer's master Xfer object so the "guido" will be set.
    int index = src->GetGuido();

    // Start at 2 because we don't want to send ViewerRPC or
    // PostponedAction back to the client.
    if(index >= 2 && index < viewerState->GetNumStateObjects())
    {
        AttributeSubject *dest = viewerState->GetStateObject(index);

        if(dest != src) // dest and src pointers are the same.
        {
            //
            // If the object is eligible for partial sends, figure out which
            // fields should be sent to the client so we can reduce the amount
            // of network traffic to the client.
            //
            intVector fieldsToSelect;
            if(viewerState->GetPartialSendFlag(index))
            {
                for(int f = 0; f < dest->NumAttributes(); ++f)
                {
                    if(src->IsSelected(f) ||
                       !dest->FieldsEqual(f, src))
                    {
                        fieldsToSelect.push_back(f);
                    }
                }
            }

            // Copy all attributes into the destination state object.
            dest->CopyAttributes(src);

            // If the object is eligible for partial sends, select the right
            // fields now.
            if(viewerState->GetPartialSendFlag(index))
                dest->SelectFields(fieldsToSelect);
        }

        // Send the state object to the client.
        SetUpdate(false);
        dest->Notify();
    }
}

//
// Qt slot functions
//

// ****************************************************************************
// Method: ViewerClientConnection::ReadFromClientAndProcess
//
// Purpose: 
//   This is a Qt slot function that reads from the client and processes
//   the results. Note that this object observes all of the state objects
//   attached to its xfer so the Update method will get called.
//
// Note:
//     This method is called as a result of the socket notifier firing.
//
// Programmer: Brad Whitlock
// Creation:   Tue May 31 13:46:34 PST 2005
//
// Modifications:
//    Brad Whitlock, Fri Jan 9 14:47:07 PST 2009
//    Added exception handling code so exceptions cannot get back into the
//    Qt event loop.
//   
// ****************************************************************************

void
ViewerClientConnection::ReadFromClientAndProcess(int)
{
    TRY
    {
        int amountRead = xfer->GetInputConnection()->Fill();

        //
        // Try and process the input.
        //
        if (amountRead > 0)
            xfer->Process();
    }
    CATCH(LostConnectionException)
    {
        // Emit a signal so the viewer can delete this dead connection.
        emit DisconnectClient(this);
    }
    CATCHALL
    {
        ; // nothing
    }
    ENDTRY
}

// ****************************************************************************
// Method: ViewerClientConnection::Name
//
// Purpose: 
//   Return the name of the client.
//
// Programmer: Brad Whitlock
// Creation:   Fri May  9 14:44:06 PDT 2008
//
// Modifications:
//   
// ****************************************************************************

const QString &
ViewerClientConnection::Name() const
{
    return name;
}

// ****************************************************************************
// Method: ViewerClientConnection::Name
//
// Purpose:
//   Return the name of the client.
//
// Programmer: Brad Whitlock
// Creation:   Fri May  9 14:44:06 PDT 2008
//
// Modifications:
//
// ****************************************************************************

void
ViewerClientConnection::ForceDisconnectClient() {
    emit DisconnectClient(this);
}







