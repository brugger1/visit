// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                    avtArrayDecomposeExpression.C                          //
// ************************************************************************* //

#include <avtArrayDecomposeExpression.h>

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cerrno>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>

#include <ExprToken.h>
#include <avtArrayMetaData.h>
#include <avtCallback.h>
#include <avtDatabase.h>
#include <avtDatabaseMetaData.h>
#include <avtExprNode.h>


#include <DebugStream.h>
#include <ExpressionException.h>
#include <ImproperUseException.h>
#include <InvalidFilesException.h>

#include <string>
#include <vector>

// ****************************************************************************
//  Method: avtArrayDecomposeExpression constructor
//
//  Purpose:
//      Defines the constructor.  Note: this should not be inlined in the
//      header because it causes problems for certain compilers.
//
//  Programmer: Hank Childs
//  Creation:   July 21, 2005
//
//  Modifications:
//
//    Alister Maguire, Tue Sep 24 11:15:10 MST 2019
//    Initialized canApplyToDirectDatabaseQOT. 
//
// ****************************************************************************

avtArrayDecomposeExpression::avtArrayDecomposeExpression()
{
    issuedWarning = false;
    canApplyToDirectDatabaseQOT = true;
    index = -1;
}


// ****************************************************************************
//  Method: avtArrayDecomposeExpression destructor
//
//  Purpose:
//      Defines the destructor.  Note: this should not be inlined in the header
//      because it causes problems for certain compilers.
//
//  Programmer: Hank Childs
//  Creation:   July 21, 2005
//
// ****************************************************************************

avtArrayDecomposeExpression::~avtArrayDecomposeExpression()
{
}


// ****************************************************************************
//  Method: avtArrayDecomposeExpression::DeriveVariable
//
//  Purpose:
//      Creates an array.
//
//  Arguments:
//      inDS      The input dataset.
//
//  Returns:      The derived variable.  The calling class must free this
//                memory.
//
//  Programmer:   Hank Childs
//  Creation:     July 21, 2005
//
//  Modifications:
//
//    Hank Childs, Fri Jun  9 14:22:43 PDT 2006
//    Remove unused variable.
//
// ****************************************************************************

vtkDataArray *
avtArrayDecomposeExpression::DeriveVariable(vtkDataSet *in_ds, int currentDomainsIndex)
{
    if (activeVariable == NULL)
        EXCEPTION2(ExpressionException, outputVariableName, 
                   "Asked to decompose an array, but did not "
                   "specify which variable to decompose");

    vtkDataArray *data = in_ds->GetPointData()->GetArray(activeVariable);
    if (data == NULL)
        data = in_ds->GetCellData()->GetArray(activeVariable);

    if (data == NULL)
        EXCEPTION2(ExpressionException, outputVariableName, 
                   "Unable to locate variable to decompose");

    if (indexStr == "")
    {
        if (index < 0)
            EXCEPTION2(ExpressionException, outputVariableName, 
                       "Index into array is not valid.");
    }

    std::string db = GetInput()->GetInfo().GetAttributes().GetFullDBName();
    ref_ptr<avtDatabase> dbp = avtCallback::GetDatabase(db, 0, NULL);
    if (*dbp == NULL)
        EXCEPTION1(InvalidFilesException, db.c_str());

    // Handle case where given index in expression may need to be mapped
    // through indices embedded in the array's meta data component names.
    // In truth, we really ought to enhance avtArrayMetaData to indicate
    // if it should be treated this way. Instead, we are using existence
    // of a compNames array of strings of the form xx...x#..# where x is
    // any non-digit and # is any digit and that the number part of the
    // names is increasing from a minimum of zero.
    size_t n;
    avtDatabaseMetaData *md = dbp->GetMetaData(currentTimeState);
    avtArrayMetaData const *amd = md->GetArray(std::string(activeVariable));
    if (amd && amd->compNames.size() &&
       (n = strcspn(amd->compNames[0].c_str(), "0123456789")) < strlen(amd->compNames[0].c_str()))
    {
        int absDistMin = INT_MAX;
        int nearestIndex = -1;
        int lastIdxVal = -1;
        bool validIdxSequence = true;
        for (size_t i = 0; i < amd->compNames.size(); i++)
        {
            // An exact match by component name wins
            if (indexStr == amd->compNames[i])
            {
                nearestIndex = i;
                break;
            }

            // get the current component index from its name and validate it
            int compIdxVal = (int) strtol(amd->compNames[i].c_str()+n, 0, 10);
            if ((errno != 0 && compIdxVal == 0) || compIdxVal < 0 || compIdxVal <= lastIdxVal)
            {
                validIdxSequence = false;
                break;
            }
            lastIdxVal = compIdxVal;

            int absDist = abs(index - compIdxVal);
            if (absDist < absDistMin)
            {
                absDistMin = absDist;
                nearestIndex = (int) i; 
            }
        }

        if (validIdxSequence && nearestIndex != -1)
            index = nearestIndex;
        else
        {
            if (indexStr != "")
            {
                EXCEPTION2(ExpressionException, outputVariableName, 
                           "Component name for array is not valid.");
            }
            else if (index >= data->GetNumberOfComponents())
            {
                EXCEPTION2(ExpressionException, outputVariableName, 
                           "Index into array is not valid.");
            }
        }
    }

    vtkDataArray *rv = data->NewInstance();
    vtkIdType nvals = data->GetNumberOfTuples();
    rv->SetNumberOfTuples(nvals);
    for (vtkIdType i = 0 ; i < nvals ; ++i)
        rv->SetTuple1(i, data->GetComponent(i, index));

    return rv;
}


// ****************************************************************************
//  Method: avtArrayDecomposeExpression::ProcessArguments
//
//  Purpose:
//      Tells the first argument to go generate itself.  Parses the second
//      argument into a list of material names.
//
//  Arguments:
//      inDS      The input dataset.
//
//  Returns:      The derived variable.  The calling class must free this
//                memory.
//
//  Programmer:   Hank Childs
//  Creation:     July 21, 2005
//
// ****************************************************************************

void
avtArrayDecomposeExpression::ProcessArguments(ArgsExpr *args, 
                                          ExprPipelineState *state)
{
    // Check the number of arguments
    std::vector<ArgExpr*> *arguments = args->GetArgs();
    size_t nargs = arguments->size();
    if (nargs != 2)
    {
        EXCEPTION2(ExpressionException, outputVariableName, 
                        "this expression must be specified with exactly two "
                        "arguments.  Usage: array_decompose(array, #)");
    }

    // Tell the first argument to create its filters.
    ArgExpr *firstarg = (*arguments)[0];
    avtExprNode *firstTree = dynamic_cast<avtExprNode*>(firstarg->GetExpr());
    firstTree->CreateFilters(state);

    ArgExpr *secondarg = (*arguments)[1];
    ExprParseTreeNode *secondTree = secondarg->GetExpr();
    std::string type = secondTree->GetTypeName();
    if (type == "IntegerConst")
        index = dynamic_cast<IntegerConstExpr*>(secondTree)->GetValue();
    else if (type == "StringConst")
        indexStr = dynamic_cast<StringConstExpr*>(secondTree)->GetValue();
    else
    {
        debug5 << "avtArrayDecomposeExpression: Second argument is not an integer index "
               << "or string component name." << endl;
        EXCEPTION2(ExpressionException, outputVariableName, "Second argument to array_decompose "
                                        "must be either a number or string.");
    }
}


// ****************************************************************************
//  Method: avtArrayDecomposeExpression::PreExecute
//
//  Purpose:
//      Called before execution.  This sets the issuedWarning flag to false.
//
//  Programmer: Hank Childs
//  Creation:   July 21, 2005
//
//  Modifications:
//    Jeremy Meredith, Thu Feb 15 12:02:51 EST 2007
//    Call inherited PreExecute before everything else.
//
// ****************************************************************************

void
avtArrayDecomposeExpression::PreExecute(void)
{
    avtSingleInputExpressionFilter::PreExecute();
    issuedWarning = false;
}


