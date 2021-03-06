// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                             avtDatasetSink.C                              //
// ************************************************************************* //

#include <avtDatasetSink.h>

#include <DebugStream.h>
#include <ImproperUseException.h>
#include <NoInputException.h>

#include <cstring>

// ****************************************************************************
//  Method: avtDatasetSink constructor
//
//  Programmer: Hank Childs
//  Creation:   May 29, 2001
//
// ****************************************************************************

avtDatasetSink::avtDatasetSink()
{
    input = NULL;
}


// ****************************************************************************
//  Method: avtDatasetSink destructor
//
//  Purpose:
//      Defines the destructor.  Note: this should not be inlined in the header
//      because it causes problems for certain compilers.
//
//  Programmer: Hank Childs
//  Creation:   February 5, 2004
//
// ****************************************************************************

avtDatasetSink::~avtDatasetSink()
{
    ;
}


// ****************************************************************************
//  Method: avtDatasetSink::SetTypedInput
//
//  Purpose:
//      Sets the input of the sink and performs some type checking.
//
//  Arguments:
//      in      The data set as a data object.
//
//  Programmer: Hank Childs
//  Creation:   May 29, 2001
//
//  Modifications:
//
//    Brad Whitlock, Thu Apr 4 14:53:33 PST 2002
//    Changed CopyTo to an inline template function.
//
//    Hank Childs, Tue Sep 10 09:02:17 PDT 2002
//    Do not assume that an input is non-NULL.
//
// ****************************************************************************

void
avtDatasetSink::SetTypedInput(avtDataObject_p in)
{
    if (*in != NULL && strcmp(in->GetType(), "avtDataset") != 0)
    {
        //
        // Should create a new exception here, but I'm under time constraints.
        //
        debug1 << "Looking for avtDataset, but found type \"" << in->GetType()
               << "\"." << endl;
        EXCEPTION0(ImproperUseException);
    }

    CopyTo(input, in);
}


// ****************************************************************************
//  Method: avtDatasetSink::GetInput
//
//  Purpose:
//      Gets the input of the sink (properly typed as an avtDataObject).
//
//  Returns:    The input of the sink.
//
//  Programmer: Hank Childs
//  Creation:   May 29, 2001
//
//  Modifications:
//    Brad Whitlock, Thu Apr 4 14:53:33 PST 2002
//    Changed CopyTo to an inline template function.
//
// ****************************************************************************

avtDataObject_p
avtDatasetSink::GetInput(void)
{
    avtDataObject_p rv;
    CopyTo(rv, input);
    return rv;
}

// ****************************************************************************
//  Method: avtDatasetSink::GetInput
//
//  Purpose:
//      Gets the input of the sink (properly typed as an avtDataObject).
//
//  Returns:    The input of the sink.
//
//  Programmer: Tom Fogal
//  Creation:   June 23, 2009
//
// ****************************************************************************
const avtDataObject_p
avtDatasetSink::GetInput(void) const
{
    avtDataObject_p rv;
    CopyTo(rv, input);
    return rv;
}


// ****************************************************************************
//  Method: avtDatasetSink::GetInputDataTree
//
//  Purpose:
//      Gets the data tree of the input.
//
//  Returns:    The avtDataTree corresponding to the input.
//
//  Programmer: Hank Childs
//  Creation:   September 19, 2000
//
// ****************************************************************************

avtDataTree_p
avtDatasetSink::GetInputDataTree()
{
    if (*input == NULL)
    {
        EXCEPTION0(NoInputException);
    }

    return input->GetDataTree();
}



