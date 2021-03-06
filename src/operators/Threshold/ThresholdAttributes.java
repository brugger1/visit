// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

package llnl.visit.operators;

import llnl.visit.ThresholdOpAttributes;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: ThresholdAttributes
//
// Purpose:
//    This class contains attributes for the threshold operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ThresholdAttributes extends ThresholdOpAttributes implements Plugin
{
    private static int ThresholdAttributes_numAdditionalAtts = 0;

    public ThresholdAttributes()
    {
        super(ThresholdAttributes_numAdditionalAtts);

    }

    public ThresholdAttributes(int nMoreFields)
    {
        super(ThresholdAttributes_numAdditionalAtts + nMoreFields);

    }

    public ThresholdAttributes(ThresholdAttributes obj)
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
        return ThresholdAttributes_numAdditionalAtts;
    }

    public boolean equals(ThresholdAttributes obj)
    {
        // Create the return value
        return (super.equals(obj) && true);
    }

    public String GetName() { return "Threshold"; }
    public String GetVersion() { return "1.0"; }

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

