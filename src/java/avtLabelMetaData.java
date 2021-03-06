// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

package llnl.visit;


// ****************************************************************************
// Class: avtLabelMetaData
//
// Purpose:
//    Contains label metadata attributes
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class avtLabelMetaData extends avtVarMetaData
{
    private static int avtLabelMetaData_numAdditionalAtts = 0;

    public avtLabelMetaData()
    {
        super(avtLabelMetaData_numAdditionalAtts);

    }

    public avtLabelMetaData(int nMoreFields)
    {
        super(avtLabelMetaData_numAdditionalAtts + nMoreFields);

    }

    public avtLabelMetaData(avtLabelMetaData obj)
    {
        super(obj);


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return avtLabelMetaData_numAdditionalAtts;
    }

    public boolean equals(avtLabelMetaData obj)
    {
        // Create the return value
        return (super.equals(obj) && true);
    }

    // Property setting methods
    // Property getting methods

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        super.WriteAtts(buf);

    }

    public void ReadAtts(int id, CommunicationBuffer buf)
    {
        super.ReadAtts(id, buf);
    }

    public String toString(String indent)
    {
        String str = new String();
        return super.toString(indent) + str;
    }


    // Attributes
}

