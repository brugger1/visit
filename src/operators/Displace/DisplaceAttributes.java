// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: DisplaceAttributes
//
// Purpose:
//    This class contains attributes for the displace operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class DisplaceAttributes extends AttributeSubject implements Plugin
{
    private static int DisplaceAttributes_numAdditionalAtts = 2;

    public DisplaceAttributes()
    {
        super(DisplaceAttributes_numAdditionalAtts);

        factor = 1;
        variable = new String("default");
    }

    public DisplaceAttributes(int nMoreFields)
    {
        super(DisplaceAttributes_numAdditionalAtts + nMoreFields);

        factor = 1;
        variable = new String("default");
    }

    public DisplaceAttributes(DisplaceAttributes obj)
    {
        super(obj);

        factor = obj.factor;
        variable = new String(obj.variable);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return DisplaceAttributes_numAdditionalAtts;
    }

    public boolean equals(DisplaceAttributes obj)
    {
        // Create the return value
        return ((factor == obj.factor) &&
                (variable.equals(obj.variable)));
    }

    public String GetName() { return "Displace"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetFactor(double factor_)
    {
        factor = factor_;
        Select(0);
    }

    public void SetVariable(String variable_)
    {
        variable = variable_;
        Select(1);
    }

    // Property getting methods
    public double GetFactor() { return factor; }
    public String GetVariable() { return variable; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDouble(factor);
        if(WriteSelect(1, buf))
            buf.WriteString(variable);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetFactor(buf.ReadDouble());
            break;
        case 1:
            SetVariable(buf.ReadString());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleToString("factor", factor, indent) + "\n";
        str = str + stringToString("variable", variable, indent) + "\n";
        return str;
    }


    // Attributes
    private double factor;
    private String variable;
}

