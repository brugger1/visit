// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                         avtLogicalAndExpression.C                         //
// ************************************************************************* //

#include <avtLogicalAndExpression.h>

#include <vtkDataArray.h>

#include <ExpressionException.h>


// ****************************************************************************
//  Method: avtLogicalAndExpression constructor
//
//  Purpose:
//      Defines the constructor.  Note: this should not be inlined in the
//      header because it causes problems for certain compilers.
//
//  Programmer: Hank Childs
//  Creation:   February 5, 2004
//
// ****************************************************************************

avtLogicalAndExpression::avtLogicalAndExpression()
{
    ;
}


// ****************************************************************************
//  Method: avtLogicalAndExpression destructor
//
//  Purpose:
//      Defines the destructor.  Note: this should not be inlined in the header
//      because it causes problems for certain compilers.
//
//  Programmer: Hank Childs
//  Creation:   February 5, 2004
//
// ****************************************************************************

avtLogicalAndExpression::~avtLogicalAndExpression()
{
    ;
}


// ****************************************************************************
//  Method: avtLogicalAndExpression::DoOperation
//
//  Purpose:
//      Takes the logical and of two arrays.
//
//  Arguments:
//      in1           The first input data array.
//      in2           The second input data array.
//      out           The output data array.
//      ncomponents   The number of components ('1' for scalar, '2' or '3' for
//                    vectors, etc.)
//      ntuples       The number of tuples (ie 'npoints' or 'ncells')
//
//  Programmer: Hank Childs
//  Creation:   August 20, 2003
//
// ****************************************************************************
 
void
avtLogicalAndExpression::DoOperation(vtkDataArray *in1, vtkDataArray *in2,
                                 vtkDataArray *out, int ncomponents,
                                 int ntuples)
{
    int in1ncomps = in1->GetNumberOfComponents();
    int in2ncomps = in2->GetNumberOfComponents();
    if (in1ncomps != 1 || in2ncomps != 1)
    {
        EXCEPTION2(ExpressionException, outputVariableName,
                   "Cannot logically and vector variables.");
    }

    for (int i = 0 ; i < ntuples ; i++)
    {
        bool val1, val2;

        if (in1->GetDataType() == VTK_UNSIGNED_CHAR)
        {
            val1 = (unsigned char) in1->GetTuple1(i);
        }
        else
        {
            val1 = (in1->GetTuple1(i) != 0. ? true : false);
        }

        if (in2->GetDataType() == VTK_UNSIGNED_CHAR)
        {
            val2 = (unsigned char) in2->GetTuple1(i);
        }
        else
        {
            val2 = (in2->GetTuple1(i) != 0. ? true : false);
        }

        unsigned char outval = val1 && val2;
        out->SetTuple1(i, outval);
    }
}


