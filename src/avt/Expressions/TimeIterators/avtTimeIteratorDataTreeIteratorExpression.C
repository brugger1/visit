// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                 avtTimeIteratorDataTreeIteratorExpression.C               //
// ************************************************************************* //

#include <avtTimeIteratorDataTreeIteratorExpression.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <DebugStream.h>
#include <ExpressionException.h>


// ****************************************************************************
//  Method: avtTimeIteratorDataTreeIteratorExpression constructor
//
//  Programmer: Hank Childs
//  Creation:   February 15, 2009
//
// ****************************************************************************

avtTimeIteratorDataTreeIteratorExpression::avtTimeIteratorDataTreeIteratorExpression()
{
    arrayIndex = 0;
}


// ****************************************************************************
//  Method: avtTimeIteratorDataTreeIteratorExpression destructor
//
//  Programmer: Hank Childs
//  Creation:   February 15, 2009
//
// ****************************************************************************

avtTimeIteratorDataTreeIteratorExpression::~avtTimeIteratorDataTreeIteratorExpression()
{
    ;
}


// ****************************************************************************
//  Method: avtTimeIteratorDataTreeIteratorExpression::ProcessDataTree
//
//  Purpose:
//      Processes a data tree.
//
//  Programmer:   Hank Childs
//  Creation:     February 15, 2009
//
// ****************************************************************************

void
avtTimeIteratorDataTreeIteratorExpression::ProcessDataTree(avtDataTree_p tree,
                                                           int ts)
{
    arrayIndex = 0;
    InternalProcessDataTree(tree, ts);
}


// ****************************************************************************
//  Method: avtTimeIteratorDataTreeIteratorExpression::InternalProcessDataTree
//
//  Purpose:
//      Processes a data tree.
//
//  Programmer:   Hank Childs
//  Creation:     February 15, 2009
//
// ****************************************************************************

void
avtTimeIteratorDataTreeIteratorExpression::InternalProcessDataTree(
                                                    avtDataTree_p tree, int ts)
{
    if (*tree == NULL)
    {
        return;
    }

    int nc = tree->GetNChildren();

    if (nc <= 0 && !tree->HasData())
    {
        return;
    }

    if (nc == 0)
    {
        //
        // there is only one dataset to process
        //
        vtkDataSet *in_ds = tree->GetDataRepresentation().GetDataVTK();
        PrepareAndExecuteDataset(in_ds, ts);
    }
    else
    {
        for (int j = 0; j < nc; j++)
            if (tree->ChildIsPresent(j))
                InternalProcessDataTree(tree->GetChild(j), ts);
    }
}


// ****************************************************************************
//  Method: avtTimeIteratorDataTreeIteratorExpression::PrepareAndExecuteDataset
//
//  Purpose:
//      Gets the proper data arrays and calls ExecuteDataset
//
//  Programmer:   Hank Childs
//  Creation:     February 15, 2009
//
// ****************************************************************************

void
avtTimeIteratorDataTreeIteratorExpression::PrepareAndExecuteDataset(
                                                        vtkDataSet *ds, int ts)
{
    std::vector<vtkDataArray *> ds_vars;
    std::vector<vtkDataArray *> delete_vars;

    bool haveZonal = false;

    int nvars = (int)varnames.size();
    if (cmfeType == POS_CMFE)
        nvars--;
    for (int i = 0 ; i < nvars ; i++)
    {
        std::string vname = GetInternalVarname(i);
        vtkDataArray *cell_data1 = ds->GetCellData()->GetArray(vname.c_str());
        vtkDataArray *point_data1 = ds->GetPointData()->GetArray(vname.c_str());
        if (cell_data1 == NULL && point_data1 == NULL)
        {
            EXCEPTION2(ExpressionException, outputVariableName,
                       "An internal error occurred when calculating an expression."
                       "  Please contact a VisIt developer.");
        }
        haveZonal = (cell_data1 != NULL);
    }

    bool doZonal = false;
    if (haveZonal)
        doZonal = true;  // mixed centering -> zonal

    for (int i = 0 ; i < nvars ; i++)
    {
        std::string vname = GetInternalVarname(i);
        vtkDataArray *cell_data1 = ds->GetCellData()->GetArray(vname.c_str());
        vtkDataArray *point_data1 = ds->GetPointData()->GetArray(vname.c_str());

        if (doZonal)
        { 
            if (cell_data1 != NULL)
                ds_vars.push_back(cell_data1);
            else
            {
                vtkDataArray *tmp = Recenter(ds, point_data1, AVT_NODECENT, 
                                             varnames[i], AVT_ZONECENT);
                ds_vars.push_back(tmp);
                delete_vars.push_back(tmp);
            }
        }
        else
        {
            if (point_data1 != NULL)
                ds_vars.push_back(point_data1);
            else
            {
                vtkDataArray *tmp = Recenter(ds, cell_data1, AVT_ZONECENT, 
                                             varnames[i], AVT_NODECENT);
                ds_vars.push_back(tmp);
                delete_vars.push_back(tmp);
            }
        }
    }

    vtkDataArray *out_arr = vars[arrayIndex++];
    ExecuteDataset(ds_vars, out_arr, ts);

    for (size_t i = 0 ; i < delete_vars.size() ; i++)
        delete_vars[i]->Delete();
}


// ****************************************************************************
//  Method: avtTimeIteratorDataTreeIteratorExpression::InitializeOutput
//
//  Purpose:
//      Sets up the data arrays we will use later.
//
//  Programmer:   Hank Childs
//  Creation:     February 15, 2009
//
//  Modifications:
//    Brad Whitlock, Mon Jul 19 14:23:18 PDT 2010
//    Clear the vars vector here too.
//
// ****************************************************************************

void
avtTimeIteratorDataTreeIteratorExpression::InitializeOutput(void)
{
    avtDataTree_p topTree = GetInputDataTree();
    std::vector<avtDataTree_p> treesToProcess;
    treesToProcess.push_back(topTree);
    size_t curIndex = 0;
    vars.clear();
    while (curIndex < treesToProcess.size())
    {
        avtDataTree_p tree = treesToProcess[curIndex];
        curIndex++;
        if (*tree == NULL)
            continue;

        int nc = tree->GetNChildren();
        if (nc <= 0 && !tree->HasData())
            continue;

        if (nc == 0)
        {
            vtkDataSet *in_ds = tree->GetDataRepresentation().GetDataVTK();
            vtkDoubleArray *arr = vtkDoubleArray::New();
            arr->SetNumberOfComponents(GetIntermediateSize());
            if (IsPointVariable())
                arr->SetNumberOfTuples(in_ds->GetNumberOfPoints());
            else
                arr->SetNumberOfTuples(in_ds->GetNumberOfCells());
            arr->SetName(outputVariableName);
            vars.push_back(arr);
        }
        else
        {
            for (int j = 0; j < nc; j++)
                if (tree->ChildIsPresent(j))
                    treesToProcess.push_back(tree->GetChild(j));
        }
    }
    
    // Set up this counter variable that we use during "ProcessDataTree".
    arrayIndex = 0;
}


// ****************************************************************************
//  Method: avtTimeIteratorDataTreeIteratorExpression::FinalizeOutput
//
//  Purpose:
//      Take the data array we have created and add it to the output.
//
//  Programmer:   Hank Childs
//  Creation:     February 15, 2009
//
// ****************************************************************************

void
avtTimeIteratorDataTreeIteratorExpression::FinalizeOutput(void)
{
    avtDataTree_p tree = GetInputDataTree();
    arrayIndex = 0;
    avtDataTree_p rv = ConstructOutput(tree);
    SetOutputDataTree(rv);
}


// ****************************************************************************
//  Method: avtTimeIteratorDataTreeIteratorExpression::ConstructOutput
//
//  Purpose:
//      A helper method for FinalizeOutput that iterates over a tree.
//
//  Programmer:   Hank Childs
//  Creation:     February 15, 2009
//
// ****************************************************************************

avtDataTree_p
avtTimeIteratorDataTreeIteratorExpression::ConstructOutput(avtDataTree_p t)
{
    if (*t == NULL)
    {
        return NULL;
    }

    int nc = t->GetNChildren();

    if (nc <= 0 && !t->HasData())
    {
        return NULL;
    }

    if (nc == 0)
    {
        //
        // there is only one dataset to process
        //
        vtkDataSet *in_ds = t->GetDataRepresentation().GetDataVTK();
        vtkDataSet *new_ds = in_ds->NewInstance();
        new_ds->ShallowCopy(in_ds);
        vtkDataArray *final_arr = 
                        ConvertIntermediateArrayToFinalArray(vars[arrayIndex]);
        vars[arrayIndex]->Delete();
        arrayIndex++;
        if (IsPointVariable())
            new_ds->GetPointData()->AddArray(final_arr);
        else
            new_ds->GetCellData()->AddArray(final_arr);
        final_arr->Delete();
        avtDataTree_p rv = new avtDataTree(new_ds,
                                   t->GetDataRepresentation().GetDomain(),
                                   t->GetDataRepresentation().GetLabel());
        new_ds->Delete();
        return rv;
    }
    else
    {
        // there is more than one input dataset to process
        // and we need an output datatree for each
        //
        avtDataTree_p *outDT = new avtDataTree_p[nc];
        for (int j = 0; j < nc; j++)
        {
            if (t->ChildIsPresent(j))
            {
                outDT[j] = ConstructOutput(t->GetChild(j));
            }
            else
            {
                outDT[j] = NULL;
            }
        }
        avtDataTree_p rv = new avtDataTree(nc, outDT);
        delete [] outDT;
        return (rv);
    }
}


// ****************************************************************************
//  avtTimeIteratorDataTreeIteratorExpression::
//                                         ConvertIntermediateArrayToFinalArray
//
//  Purpose:
//      Some derived types may use a data array when doing its processing that
//      is bigger than the output.  If that is the case, we need to process
//      the output before exiting.  It is up to the derived type to do this.
//      So this method is a hook ... basically a no-op ... for those derived
//      types.
//
//  Note:   It is assumed that the return value of this method should have its
//          reference count decremented.
//
//  Programmer: Hank Childs
//  Creation:   February 15, 2009
//
// ****************************************************************************

vtkDataArray *
avtTimeIteratorDataTreeIteratorExpression::ConvertIntermediateArrayToFinalArray
                                             (vtkDataArray *intermediateArray)
{
    if (GetIntermediateSize() != GetVariableDimension())
    {
        // If the intermediate size is redefined, then this method must be
        // redefined as well.
        EXCEPTION0(ImproperUseException);
    }

    intermediateArray->Register(NULL);
    return intermediateArray;
}


// ****************************************************************************
//  Method:  avtTimeIteratorDataTreeIteratorExpression::GetVariableType
//
//  Purpose:
//    Try to do better than unknown type for the output type.
//    Specifically, the time iteration stuff typically makes the output
//    num components the same as the input, so we try to re-use the
//    input data type.
//
//  Arguments:
//    
//
//  Programmer:  Jeremy Meredith
//  Creation:    March 18, 2009
//
// ****************************************************************************
avtVarType
avtTimeIteratorDataTreeIteratorExpression::GetVariableType()
{
    if (varnames.size() != 1)
        return AVT_UNKNOWN_TYPE;

    avtDataAttributes &inatts = GetInput()->GetInfo().GetAttributes();
    if (!inatts.ValidVariable(varnames[0]))
        return AVT_UNKNOWN_TYPE;

    return inatts.GetVariableType(varnames[0].c_str());
}

// ****************************************************************************
//  Method:  avtTimeIteratorDataTreeIteratorExpression::GetVariableType
//
//  Purpose:
//    Try to do better than use the "primary" variable for the dimension
//    since we're not getting the right primary variable.  Similar to 
//    GetVariableType, our output dimension should be the same as the
//    input dimension.
//
//  Arguments:
//    
//
//  Programmer:  Jeremy Meredith
//  Creation:    March 18, 2009
//
// ****************************************************************************
int
avtTimeIteratorDataTreeIteratorExpression::GetVariableDimension(void)
{
    if (varnames.size() != 1)
        return 1;

    avtDataAttributes &inatts = GetInput()->GetInfo().GetAttributes();
    if (!inatts.ValidVariable(varnames[0]))
        return 1;

    return inatts.GetVariableDimension(varnames[0].c_str());
}

