// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                             avtMinMaxExpression.C                         //
// ************************************************************************* //

#include <avtMinMaxExpression.h>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>

#include <ExpressionException.h>
#include <DebugStream.h>


// ****************************************************************************
//  Method: avtMinMaxExpression constructor
//
//  Purpose:
//      Defines the constructor.  Note: this should not be inlined in the
//      header because it causes problems for certain compilers.
//
//  Programmer: Hank Childs
//  Creation:   March 13, 2006
//
//  Modifications:
//
//      Eddie Rusu, Mon Sep 30 14:49:38 PDT 2019
//      Added alternate constructor that takes doMin as input.
//
// ****************************************************************************

avtMinMaxExpression::avtMinMaxExpression()
{
    doMin = false;
}

avtMinMaxExpression::avtMinMaxExpression(bool _doMin)
{
    doMin = _doMin;
}


// ****************************************************************************
//  Method: avtMinMaxExpression destructor
//
//  Purpose:
//      Defines the destructor.  Note: this should not be inlined in the header
//      because it causes problems for certain compilers.
//
//  Programmer: Hank Childs
//  Creation:   March 13, 2006
//
// ****************************************************************************

avtMinMaxExpression::~avtMinMaxExpression()
{
    ;
}


// ****************************************************************************
//  Method: avtMinMaxExpression::DoOperation
//
//  Purpose:
//      Finds the max or min among all the input variables stored in
//      dataArrays.
//
//  Programmer: Eddie Rusu
//  Creation:   Mon Sep 30 14:09:49 PDT 2019
//
// ****************************************************************************

vtkDataArray*
avtMinMaxExpression::DoOperation()
{
    debug4 << "Entering avtMinMaxExpression::DoOperation()" << std::endl;
    vtkDataArray* output = CreateOutputVariable();

    // Loop over all inputs and determine the min/max
    DoOperationHelper(output, dataArrays[0], dataArrays[1]);
    for (int i = 2; i < nProcessedArgs; ++i)
    {
        DoOperationHelper(output, output, dataArrays[i]);
    }

    debug4 << "Exiting  avtMinMaxExpression::DoOperation()" << std::endl;
    return output;
}

// ****************************************************************************
//  Method: avtMinMaxExpression::DoOperationHelper
//
//  Purpose:
//      Finds the minimum or maximum value between two arrays
//
//  Arguments:
//      out     The output variable.
//      in1     The first input data array.
//      in2     The second input data array.
//
//  Programmer: Hank Childs
//  Creation:   March 13, 2006
//
//  Modifications:
//
//    Hank Childs, Mon Jan 14 18:26:58 PST 2008
//    Add support for singleton constants.
//
//    Kathleen Biagas, Wed Apr 4 12:13:10 PDT 2012
//    Change float to double.
//
//    Eddie Rusu, Mon Sep 30 11:13:18 PDT 2019
//    Modified inputs for use with new multi-variable architecture. Changed
//    name from DoOperation to DoOperationHelper.
//
// ****************************************************************************

void
avtMinMaxExpression::DoOperationHelper(vtkDataArray *out, vtkDataArray *in1,
        vtkDataArray *in2)
{
    debug5 << "Entering avtMinMaxExpression::DoOperationHelper()" << std::endl;
    bool var1IsSingleton = (in1->GetNumberOfTuples() == 1);
    bool var2IsSingleton = (in2->GetNumberOfTuples() == 1);

    int ntuples = out->GetNumberOfTuples();
    int in1ncomps = in1->GetNumberOfComponents();
    int in2ncomps = in2->GetNumberOfComponents();

    if (in1ncomps == in2ncomps)
    {
        for (int i = 0 ; i < ntuples ; i++)
        {
            for (int j = 0 ; j < in1ncomps ; j++)
            {
                vtkIdType tup1 = (var1IsSingleton ? 0 : i);
                vtkIdType tup2 = (var2IsSingleton ? 0 : i);
                double val1 = in1->GetComponent(tup1, j);
                double val2 = in2->GetComponent(tup2, j);
                bool val1Bigger = (val1 > val2);
                // Circumflex (^) is the exclusive or.
                // doMin == true  && val1Bigger == true  --> val2
                // doMin == false && val1Bigger == true  --> val1
                // doMin == true  && val1Bigger == false --> val1
                // doMin == false && val1Bigger == false --> val2
                //  --> values same, then val2, values different, then val1
                double outval = (doMin ^ val1Bigger ? val1 : val2);
                out->SetComponent(i, j, outval);
            }
        }
    }
    else if (in1ncomps > 1 && in2ncomps == 1)
    {
        for (int i = 0 ; i < ntuples ; i++)
        {
            vtkIdType tup1 = (var1IsSingleton ? 0 : i);
            vtkIdType tup2 = (var2IsSingleton ? 0 : i);
            double val2 = in2->GetTuple1(tup2);
            for (int j = 0 ; j < in1ncomps ; j++)
            {
                double val1 = in1->GetComponent(tup1, j);
                bool val1Bigger = (val1 > val2);
                // Circumflex (^) is the exclusive or.
                // doMin == true  && val1Bigger == true  --> val2
                // doMin == false && val1Bigger == true  --> val1
                // doMin == true  && val1Bigger == false --> val1
                // doMin == false && val1Bigger == false --> val2
                //  --> values same, then val2, values different, then val1
                double outval = (doMin ^ val1Bigger ? val1 : val2);
                out->SetComponent(i, j, outval);
            }
        }
    }
    else if (in1ncomps == 1 && in2ncomps > 1)
    {
        for (int i = 0 ; i < ntuples ; i++)
        {
            vtkIdType tup1 = (var1IsSingleton ? 0 : i);
            vtkIdType tup2 = (var2IsSingleton ? 0 : i);
            double val1 = in1->GetTuple1(tup1);
            for (int j = 0 ; j < in2ncomps ; j++)
            {
                double val2 = in2->GetComponent(tup2, j);
                bool val1Bigger = (val1 > val2);
                // Circumflex (^) is the exclusive or.
                // doMin == true  && val1Bigger == true  --> val2
                // doMin == false && val1Bigger == true  --> val1
                // doMin == true  && val1Bigger == false --> val1
                // doMin == false && val1Bigger == false --> val2
                //  --> values same, then val2, values different, then val1
                double outval = (doMin ^ val1Bigger ? val1 : val2);
                out->SetComponent(i, j, outval);
            }
        }
    }
    else
    {
        EXCEPTION2(ExpressionException, outputVariableName, 
                "Don't know how to take minimums or maximums with data of "
                "differing dimensions.");
    }
    debug5 << "Exiting  avtMinMaxExpression::DoOperationHelper()" << std::endl;
}
