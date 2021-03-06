// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

package llnl.visit;

import java.util.Vector;

// ****************************************************************************
// Class: ClientInformationList
//
// Purpose:
//    Contains the information for all connected clients.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ClientInformationList extends AttributeSubject
{
    private static int ClientInformationList_numAdditionalAtts = 1;

    public ClientInformationList()
    {
        super(ClientInformationList_numAdditionalAtts);

        clients = new Vector();
    }

    public ClientInformationList(int nMoreFields)
    {
        super(ClientInformationList_numAdditionalAtts + nMoreFields);

        clients = new Vector();
    }

    public ClientInformationList(ClientInformationList obj)
    {
        super(obj);

        int i;

        // *** Copy the clients field ***
        clients = new Vector(obj.clients.size());
        for(i = 0; i < obj.clients.size(); ++i)
        {
            ClientInformation oldObj = (ClientInformation)obj.clients.elementAt(i);
            clients.addElement(new ClientInformation(oldObj));
        }


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ClientInformationList_numAdditionalAtts;
    }

    public boolean equals(ClientInformationList obj)
    {
        int i;

        // Compare the elements in the clients vector.
        boolean clients_equal = (obj.clients.size() == clients.size());
        for(i = 0; (i < clients.size()) && clients_equal; ++i)
        {
            // Make references to ClientInformation from Object.
            ClientInformation clients1 = (ClientInformation)clients.elementAt(i);
            ClientInformation clients2 = (ClientInformation)obj.clients.elementAt(i);
            clients_equal = clients1.equals(clients2);
        }
        // Create the return value
        return (clients_equal);
    }

    // Property setting methods
    // Property getting methods
    public Vector GetClients() { return clients; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
        {
            buf.WriteInt(clients.size());
            for(int i = 0; i < clients.size(); ++i)
            {
                ClientInformation tmp = (ClientInformation)clients.elementAt(i);
                tmp.Write(buf);
            }
        }
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        {
            int len = buf.ReadInt();
            clients.clear();
            for(int j = 0; j < len; ++j)
            {
                ClientInformation tmp = new ClientInformation();
                tmp.Read(buf);
                clients.addElement(tmp);
            }
        }
        Select(0);
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "clients = {\n";
        for(int i = 0; i < clients.size(); ++i)
        {
            AttributeSubject s = (AttributeSubject)clients.elementAt(i);
            str = str + s.toString(indent + "    ");
            if(i < clients.size()-1)
                str = str + ", ";
            str = str + "\n";
        }
        str = str + "}\n";
        return str;
    }

    // Attributegroup convenience methods
    public void AddClients(ClientInformation obj)
    {
        clients.addElement(new ClientInformation(obj));
        Select(0);
    }

    public void ClearClients()
    {
        clients.clear();
        Select(0);
    }

    public void RemoveClients(int index)
    {
        if(index >= 0 && index < clients.size())
        {
            clients.remove(index);
            Select(0);
        }
    }

    public int GetNumClients()
    {
        return clients.size();
    }

    public ClientInformation GetClients(int i)
    {
        ClientInformation tmp = (ClientInformation)clients.elementAt(i);
        return tmp;
    }


    // Attributes
    private Vector clients; // vector of ClientInformation objects
}

