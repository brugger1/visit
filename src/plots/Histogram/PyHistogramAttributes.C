// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include <PyHistogramAttributes.h>
#include <ObserverToCallback.h>
#include <stdio.h>
#include <ColorAttribute.h>

// ****************************************************************************
// Module: PyHistogramAttributes
//
// Purpose: 
//   Attributes for Histogram Plot
//
// Note:       Autogenerated by xml2python. Do not modify by hand!
//
// Programmer: xml2python
// Creation:   omitted
//
// ****************************************************************************

//
// This struct contains the Python type information and a HistogramAttributes.
//
struct HistogramAttributesObject
{
    PyObject_HEAD
    HistogramAttributes *data;
    bool        owns;
    PyObject   *parent;
};

//
// Internal prototypes
//
static PyObject *NewHistogramAttributes(int);

std::string
PyHistogramAttributes_ToString(const HistogramAttributes *atts, const char *prefix)
{
    std::string str;
    char tmpStr[1000];

    const char *basedOn_names = "ManyVarsForSingleZone, ManyZonesForSingleVar";
    switch (atts->GetBasedOn())
    {
      case HistogramAttributes::ManyVarsForSingleZone:
          snprintf(tmpStr, 1000, "%sbasedOn = %sManyVarsForSingleZone  # %s\n", prefix, prefix, basedOn_names);
          str += tmpStr;
          break;
      case HistogramAttributes::ManyZonesForSingleVar:
          snprintf(tmpStr, 1000, "%sbasedOn = %sManyZonesForSingleVar  # %s\n", prefix, prefix, basedOn_names);
          str += tmpStr;
          break;
      default:
          break;
    }

    const char *histogramType_names = "Frequency, Weighted, Variable";
    switch (atts->GetHistogramType())
    {
      case HistogramAttributes::Frequency:
          snprintf(tmpStr, 1000, "%shistogramType = %sFrequency  # %s\n", prefix, prefix, histogramType_names);
          str += tmpStr;
          break;
      case HistogramAttributes::Weighted:
          snprintf(tmpStr, 1000, "%shistogramType = %sWeighted  # %s\n", prefix, prefix, histogramType_names);
          str += tmpStr;
          break;
      case HistogramAttributes::Variable:
          snprintf(tmpStr, 1000, "%shistogramType = %sVariable  # %s\n", prefix, prefix, histogramType_names);
          str += tmpStr;
          break;
      default:
          break;
    }

    snprintf(tmpStr, 1000, "%sweightVariable = \"%s\"\n", prefix, atts->GetWeightVariable().c_str());
    str += tmpStr;
    const char *limitsMode_names = "OriginalData, CurrentPlot";
    switch (atts->GetLimitsMode())
    {
      case HistogramAttributes::OriginalData:
          snprintf(tmpStr, 1000, "%slimitsMode = %sOriginalData  # %s\n", prefix, prefix, limitsMode_names);
          str += tmpStr;
          break;
      case HistogramAttributes::CurrentPlot:
          snprintf(tmpStr, 1000, "%slimitsMode = %sCurrentPlot  # %s\n", prefix, prefix, limitsMode_names);
          str += tmpStr;
          break;
      default:
          break;
    }

    if(atts->GetMinFlag())
        snprintf(tmpStr, 1000, "%sminFlag = 1\n", prefix);
    else
        snprintf(tmpStr, 1000, "%sminFlag = 0\n", prefix);
    str += tmpStr;
    if(atts->GetMaxFlag())
        snprintf(tmpStr, 1000, "%smaxFlag = 1\n", prefix);
    else
        snprintf(tmpStr, 1000, "%smaxFlag = 0\n", prefix);
    str += tmpStr;
    snprintf(tmpStr, 1000, "%smin = %g\n", prefix, atts->GetMin());
    str += tmpStr;
    snprintf(tmpStr, 1000, "%smax = %g\n", prefix, atts->GetMax());
    str += tmpStr;
    snprintf(tmpStr, 1000, "%snumBins = %d\n", prefix, atts->GetNumBins());
    str += tmpStr;
    snprintf(tmpStr, 1000, "%sdomain = %d\n", prefix, atts->GetDomain());
    str += tmpStr;
    snprintf(tmpStr, 1000, "%szone = %d\n", prefix, atts->GetZone());
    str += tmpStr;
    if(atts->GetUseBinWidths())
        snprintf(tmpStr, 1000, "%suseBinWidths = 1\n", prefix);
    else
        snprintf(tmpStr, 1000, "%suseBinWidths = 0\n", prefix);
    str += tmpStr;
    const char *outputType_names = "Curve, Block";
    switch (atts->GetOutputType())
    {
      case HistogramAttributes::Curve:
          snprintf(tmpStr, 1000, "%soutputType = %sCurve  # %s\n", prefix, prefix, outputType_names);
          str += tmpStr;
          break;
      case HistogramAttributes::Block:
          snprintf(tmpStr, 1000, "%soutputType = %sBlock  # %s\n", prefix, prefix, outputType_names);
          str += tmpStr;
          break;
      default:
          break;
    }

    snprintf(tmpStr, 1000, "%slineWidth = %d\n", prefix, atts->GetLineWidth());
    str += tmpStr;
    const unsigned char *color = atts->GetColor().GetColor();
    snprintf(tmpStr, 1000, "%scolor = (%d, %d, %d, %d)\n", prefix, int(color[0]), int(color[1]), int(color[2]), int(color[3]));
    str += tmpStr;
    const char *dataScale_names = "Linear, Log, SquareRoot";
    switch (atts->GetDataScale())
    {
      case HistogramAttributes::Linear:
          snprintf(tmpStr, 1000, "%sdataScale = %sLinear  # %s\n", prefix, prefix, dataScale_names);
          str += tmpStr;
          break;
      case HistogramAttributes::Log:
          snprintf(tmpStr, 1000, "%sdataScale = %sLog  # %s\n", prefix, prefix, dataScale_names);
          str += tmpStr;
          break;
      case HistogramAttributes::SquareRoot:
          snprintf(tmpStr, 1000, "%sdataScale = %sSquareRoot  # %s\n", prefix, prefix, dataScale_names);
          str += tmpStr;
          break;
      default:
          break;
    }

    const char *binScale_names = "Linear, Log, SquareRoot";
    switch (atts->GetBinScale())
    {
      case HistogramAttributes::Linear:
          snprintf(tmpStr, 1000, "%sbinScale = %sLinear  # %s\n", prefix, prefix, binScale_names);
          str += tmpStr;
          break;
      case HistogramAttributes::Log:
          snprintf(tmpStr, 1000, "%sbinScale = %sLog  # %s\n", prefix, prefix, binScale_names);
          str += tmpStr;
          break;
      case HistogramAttributes::SquareRoot:
          snprintf(tmpStr, 1000, "%sbinScale = %sSquareRoot  # %s\n", prefix, prefix, binScale_names);
          str += tmpStr;
          break;
      default:
          break;
    }

    if(atts->GetNormalizeHistogram())
        snprintf(tmpStr, 1000, "%snormalizeHistogram = 1\n", prefix);
    else
        snprintf(tmpStr, 1000, "%snormalizeHistogram = 0\n", prefix);
    str += tmpStr;
    if(atts->GetComputeAsCDF())
        snprintf(tmpStr, 1000, "%scomputeAsCDF = 1\n", prefix);
    else
        snprintf(tmpStr, 1000, "%scomputeAsCDF = 0\n", prefix);
    str += tmpStr;
    return str;
}

static PyObject *
HistogramAttributes_Notify(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    obj->data->Notify();
    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_SetBasedOn(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the basedOn in the object.
    if(ival >= 0 && ival < 2)
        obj->data->SetBasedOn(HistogramAttributes::BasedOn(ival));
    else
    {
        fprintf(stderr, "An invalid basedOn value was given. "
                        "Valid values are in the range of [0,1]. "
                        "You can also use the following names: "
                        "ManyVarsForSingleZone, ManyZonesForSingleVar.");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetBasedOn(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(long(obj->data->GetBasedOn()));
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetHistogramType(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the histogramType in the object.
    if(ival >= 0 && ival < 3)
        obj->data->SetHistogramType(HistogramAttributes::BinContribution(ival));
    else
    {
        fprintf(stderr, "An invalid histogramType value was given. "
                        "Valid values are in the range of [0,2]. "
                        "You can also use the following names: "
                        "Frequency, Weighted, Variable.");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetHistogramType(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(long(obj->data->GetHistogramType()));
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetWeightVariable(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    char *str;
    if(!PyArg_ParseTuple(args, "s", &str))
        return NULL;

    // Set the weightVariable in the object.
    obj->data->SetWeightVariable(std::string(str));

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetWeightVariable(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyString_FromString(obj->data->GetWeightVariable().c_str());
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetLimitsMode(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the limitsMode in the object.
    if(ival >= 0 && ival < 2)
        obj->data->SetLimitsMode(HistogramAttributes::LimitsMode(ival));
    else
    {
        fprintf(stderr, "An invalid limitsMode value was given. "
                        "Valid values are in the range of [0,1]. "
                        "You can also use the following names: "
                        "OriginalData, CurrentPlot.");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetLimitsMode(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(long(obj->data->GetLimitsMode()));
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetMinFlag(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the minFlag in the object.
    obj->data->SetMinFlag(ival != 0);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetMinFlag(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(obj->data->GetMinFlag()?1L:0L);
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetMaxFlag(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the maxFlag in the object.
    obj->data->SetMaxFlag(ival != 0);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetMaxFlag(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(obj->data->GetMaxFlag()?1L:0L);
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetMin(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    double dval;
    if(!PyArg_ParseTuple(args, "d", &dval))
        return NULL;

    // Set the min in the object.
    obj->data->SetMin(dval);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetMin(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyFloat_FromDouble(obj->data->GetMin());
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetMax(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    double dval;
    if(!PyArg_ParseTuple(args, "d", &dval))
        return NULL;

    // Set the max in the object.
    obj->data->SetMax(dval);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetMax(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyFloat_FromDouble(obj->data->GetMax());
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetNumBins(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the numBins in the object.
    obj->data->SetNumBins((int)ival);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetNumBins(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(long(obj->data->GetNumBins()));
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetDomain(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the domain in the object.
    obj->data->SetDomain((int)ival);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetDomain(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(long(obj->data->GetDomain()));
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetZone(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the zone in the object.
    obj->data->SetZone((int)ival);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetZone(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(long(obj->data->GetZone()));
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetUseBinWidths(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the useBinWidths in the object.
    obj->data->SetUseBinWidths(ival != 0);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetUseBinWidths(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(obj->data->GetUseBinWidths()?1L:0L);
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetOutputType(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the outputType in the object.
    if(ival >= 0 && ival < 2)
        obj->data->SetOutputType(HistogramAttributes::OutputType(ival));
    else
    {
        fprintf(stderr, "An invalid outputType value was given. "
                        "Valid values are in the range of [0,1]. "
                        "You can also use the following names: "
                        "Curve, Block.");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetOutputType(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(long(obj->data->GetOutputType()));
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetLineWidth(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the lineWidth in the object.
    obj->data->SetLineWidth(ival);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetLineWidth(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(long(obj->data->GetLineWidth()));
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetColor(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int c[4];
    if(!PyArg_ParseTuple(args, "iiii", &c[0], &c[1], &c[2], &c[3]))
    {
        c[3] = 255;
        if(!PyArg_ParseTuple(args, "iii", &c[0], &c[1], &c[2]))
        {
            double dr, dg, db, da;
            if(PyArg_ParseTuple(args, "dddd", &dr, &dg, &db, &da))
            {
                c[0] = int(dr);
                c[1] = int(dg);
                c[2] = int(db);
                c[3] = int(da);
            }
            else if(PyArg_ParseTuple(args, "ddd", &dr, &dg, &db))
            {
                c[0] = int(dr);
                c[1] = int(dg);
                c[2] = int(db);
                c[3] = 255;
            }
            else
            {
                PyObject *tuple = NULL;
                if(!PyArg_ParseTuple(args, "O", &tuple))
                    return NULL;

                if(!PyTuple_Check(tuple))
                    return NULL;

                // Make sure that the tuple is the right size.
                if(PyTuple_Size(tuple) < 3 || PyTuple_Size(tuple) > 4)
                    return NULL;

                // Make sure that all elements in the tuple are ints.
                for(int i = 0; i < PyTuple_Size(tuple); ++i)
                {
                    PyObject *item = PyTuple_GET_ITEM(tuple, i);
                    if(PyInt_Check(item))
                        c[i] = int(PyInt_AS_LONG(PyTuple_GET_ITEM(tuple, i)));
                    else if(PyFloat_Check(item))
                        c[i] = int(PyFloat_AS_DOUBLE(PyTuple_GET_ITEM(tuple, i)));
                    else
                        return NULL;
                }
            }
        }
        PyErr_Clear();
    }

    // Set the color in the object.
    ColorAttribute ca(c[0], c[1], c[2], c[3]);
    obj->data->SetColor(ca);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetColor(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    // Allocate a tuple the with enough entries to hold the color.
    PyObject *retval = PyTuple_New(4);
    const unsigned char *color = obj->data->GetColor().GetColor();
    PyTuple_SET_ITEM(retval, 0, PyInt_FromLong(long(color[0])));
    PyTuple_SET_ITEM(retval, 1, PyInt_FromLong(long(color[1])));
    PyTuple_SET_ITEM(retval, 2, PyInt_FromLong(long(color[2])));
    PyTuple_SET_ITEM(retval, 3, PyInt_FromLong(long(color[3])));
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetDataScale(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the dataScale in the object.
    if(ival >= 0 && ival < 3)
        obj->data->SetDataScale(HistogramAttributes::DataScale(ival));
    else
    {
        fprintf(stderr, "An invalid dataScale value was given. "
                        "Valid values are in the range of [0,2]. "
                        "You can also use the following names: "
                        "Linear, Log, SquareRoot.");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetDataScale(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(long(obj->data->GetDataScale()));
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetBinScale(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the binScale in the object.
    if(ival >= 0 && ival < 3)
        obj->data->SetBinScale(HistogramAttributes::DataScale(ival));
    else
    {
        fprintf(stderr, "An invalid binScale value was given. "
                        "Valid values are in the range of [0,2]. "
                        "You can also use the following names: "
                        "Linear, Log, SquareRoot.");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetBinScale(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(long(obj->data->GetBinScale()));
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetNormalizeHistogram(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the normalizeHistogram in the object.
    obj->data->SetNormalizeHistogram(ival != 0);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetNormalizeHistogram(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(obj->data->GetNormalizeHistogram()?1L:0L);
    return retval;
}

/*static*/ PyObject *
HistogramAttributes_SetComputeAsCDF(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;

    int ival;
    if(!PyArg_ParseTuple(args, "i", &ival))
        return NULL;

    // Set the computeAsCDF in the object.
    obj->data->SetComputeAsCDF(ival != 0);

    Py_INCREF(Py_None);
    return Py_None;
}

/*static*/ PyObject *
HistogramAttributes_GetComputeAsCDF(PyObject *self, PyObject *args)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)self;
    PyObject *retval = PyInt_FromLong(obj->data->GetComputeAsCDF()?1L:0L);
    return retval;
}



PyMethodDef PyHistogramAttributes_methods[HISTOGRAMATTRIBUTES_NMETH] = {
    {"Notify", HistogramAttributes_Notify, METH_VARARGS},
    {"SetBasedOn", HistogramAttributes_SetBasedOn, METH_VARARGS},
    {"GetBasedOn", HistogramAttributes_GetBasedOn, METH_VARARGS},
    {"SetHistogramType", HistogramAttributes_SetHistogramType, METH_VARARGS},
    {"GetHistogramType", HistogramAttributes_GetHistogramType, METH_VARARGS},
    {"SetWeightVariable", HistogramAttributes_SetWeightVariable, METH_VARARGS},
    {"GetWeightVariable", HistogramAttributes_GetWeightVariable, METH_VARARGS},
    {"SetLimitsMode", HistogramAttributes_SetLimitsMode, METH_VARARGS},
    {"GetLimitsMode", HistogramAttributes_GetLimitsMode, METH_VARARGS},
    {"SetMinFlag", HistogramAttributes_SetMinFlag, METH_VARARGS},
    {"GetMinFlag", HistogramAttributes_GetMinFlag, METH_VARARGS},
    {"SetMaxFlag", HistogramAttributes_SetMaxFlag, METH_VARARGS},
    {"GetMaxFlag", HistogramAttributes_GetMaxFlag, METH_VARARGS},
    {"SetMin", HistogramAttributes_SetMin, METH_VARARGS},
    {"GetMin", HistogramAttributes_GetMin, METH_VARARGS},
    {"SetMax", HistogramAttributes_SetMax, METH_VARARGS},
    {"GetMax", HistogramAttributes_GetMax, METH_VARARGS},
    {"SetNumBins", HistogramAttributes_SetNumBins, METH_VARARGS},
    {"GetNumBins", HistogramAttributes_GetNumBins, METH_VARARGS},
    {"SetDomain", HistogramAttributes_SetDomain, METH_VARARGS},
    {"GetDomain", HistogramAttributes_GetDomain, METH_VARARGS},
    {"SetZone", HistogramAttributes_SetZone, METH_VARARGS},
    {"GetZone", HistogramAttributes_GetZone, METH_VARARGS},
    {"SetUseBinWidths", HistogramAttributes_SetUseBinWidths, METH_VARARGS},
    {"GetUseBinWidths", HistogramAttributes_GetUseBinWidths, METH_VARARGS},
    {"SetOutputType", HistogramAttributes_SetOutputType, METH_VARARGS},
    {"GetOutputType", HistogramAttributes_GetOutputType, METH_VARARGS},
    {"SetLineWidth", HistogramAttributes_SetLineWidth, METH_VARARGS},
    {"GetLineWidth", HistogramAttributes_GetLineWidth, METH_VARARGS},
    {"SetColor", HistogramAttributes_SetColor, METH_VARARGS},
    {"GetColor", HistogramAttributes_GetColor, METH_VARARGS},
    {"SetDataScale", HistogramAttributes_SetDataScale, METH_VARARGS},
    {"GetDataScale", HistogramAttributes_GetDataScale, METH_VARARGS},
    {"SetBinScale", HistogramAttributes_SetBinScale, METH_VARARGS},
    {"GetBinScale", HistogramAttributes_GetBinScale, METH_VARARGS},
    {"SetNormalizeHistogram", HistogramAttributes_SetNormalizeHistogram, METH_VARARGS},
    {"GetNormalizeHistogram", HistogramAttributes_GetNormalizeHistogram, METH_VARARGS},
    {"SetComputeAsCDF", HistogramAttributes_SetComputeAsCDF, METH_VARARGS},
    {"GetComputeAsCDF", HistogramAttributes_GetComputeAsCDF, METH_VARARGS},
    {NULL, NULL}
};

//
// Type functions
//

static void
HistogramAttributes_dealloc(PyObject *v)
{
   HistogramAttributesObject *obj = (HistogramAttributesObject *)v;
   if(obj->parent != 0)
       Py_DECREF(obj->parent);
   if(obj->owns)
       delete obj->data;
}

static int
HistogramAttributes_compare(PyObject *v, PyObject *w)
{
    HistogramAttributes *a = ((HistogramAttributesObject *)v)->data;
    HistogramAttributes *b = ((HistogramAttributesObject *)w)->data;
    return (*a == *b) ? 0 : -1;
}

PyObject *
PyHistogramAttributes_getattr(PyObject *self, char *name)
{
    if(strcmp(name, "basedOn") == 0)
        return HistogramAttributes_GetBasedOn(self, NULL);
    if(strcmp(name, "ManyVarsForSingleZone") == 0)
        return PyInt_FromLong(long(HistogramAttributes::ManyVarsForSingleZone));
    if(strcmp(name, "ManyZonesForSingleVar") == 0)
        return PyInt_FromLong(long(HistogramAttributes::ManyZonesForSingleVar));

    if(strcmp(name, "histogramType") == 0)
        return HistogramAttributes_GetHistogramType(self, NULL);
    if(strcmp(name, "Frequency") == 0)
        return PyInt_FromLong(long(HistogramAttributes::Frequency));
    if(strcmp(name, "Weighted") == 0)
        return PyInt_FromLong(long(HistogramAttributes::Weighted));
    if(strcmp(name, "Variable") == 0)
        return PyInt_FromLong(long(HistogramAttributes::Variable));

    if(strcmp(name, "weightVariable") == 0)
        return HistogramAttributes_GetWeightVariable(self, NULL);
    if(strcmp(name, "limitsMode") == 0)
        return HistogramAttributes_GetLimitsMode(self, NULL);
    if(strcmp(name, "OriginalData") == 0)
        return PyInt_FromLong(long(HistogramAttributes::OriginalData));
    if(strcmp(name, "CurrentPlot") == 0)
        return PyInt_FromLong(long(HistogramAttributes::CurrentPlot));

    if(strcmp(name, "minFlag") == 0)
        return HistogramAttributes_GetMinFlag(self, NULL);
    if(strcmp(name, "maxFlag") == 0)
        return HistogramAttributes_GetMaxFlag(self, NULL);
    if(strcmp(name, "min") == 0)
        return HistogramAttributes_GetMin(self, NULL);
    if(strcmp(name, "max") == 0)
        return HistogramAttributes_GetMax(self, NULL);
    if(strcmp(name, "numBins") == 0)
        return HistogramAttributes_GetNumBins(self, NULL);
    if(strcmp(name, "domain") == 0)
        return HistogramAttributes_GetDomain(self, NULL);
    if(strcmp(name, "zone") == 0)
        return HistogramAttributes_GetZone(self, NULL);
    if(strcmp(name, "useBinWidths") == 0)
        return HistogramAttributes_GetUseBinWidths(self, NULL);
    if(strcmp(name, "outputType") == 0)
        return HistogramAttributes_GetOutputType(self, NULL);
    if(strcmp(name, "Curve") == 0)
        return PyInt_FromLong(long(HistogramAttributes::Curve));
    if(strcmp(name, "Block") == 0)
        return PyInt_FromLong(long(HistogramAttributes::Block));

    if(strcmp(name, "lineWidth") == 0)
        return HistogramAttributes_GetLineWidth(self, NULL);
    if(strcmp(name, "color") == 0)
        return HistogramAttributes_GetColor(self, NULL);
    if(strcmp(name, "dataScale") == 0)
        return HistogramAttributes_GetDataScale(self, NULL);
    if(strcmp(name, "Linear") == 0)
        return PyInt_FromLong(long(HistogramAttributes::Linear));
    if(strcmp(name, "Log") == 0)
        return PyInt_FromLong(long(HistogramAttributes::Log));
    if(strcmp(name, "SquareRoot") == 0)
        return PyInt_FromLong(long(HistogramAttributes::SquareRoot));

    if(strcmp(name, "binScale") == 0)
        return HistogramAttributes_GetBinScale(self, NULL);
    if(strcmp(name, "Linear") == 0)
        return PyInt_FromLong(long(HistogramAttributes::Linear));
    if(strcmp(name, "Log") == 0)
        return PyInt_FromLong(long(HistogramAttributes::Log));
    if(strcmp(name, "SquareRoot") == 0)
        return PyInt_FromLong(long(HistogramAttributes::SquareRoot));

    if(strcmp(name, "normalizeHistogram") == 0)
        return HistogramAttributes_GetNormalizeHistogram(self, NULL);
    if(strcmp(name, "computeAsCDF") == 0)
        return HistogramAttributes_GetComputeAsCDF(self, NULL);

    // Try and handle legacy fields

    // lineStyle and it's possible enumerations
    bool lineStyleFound = false;
    if (strcmp(name, "lineStyle") == 0)
    {
        lineStyleFound = true;
    }
    else if (strcmp(name, "SOLID") == 0)
    {
        lineStyleFound = true;
    }
    else if (strcmp(name, "DASH") == 0)
    {
        lineStyleFound = true;
    }
    else if (strcmp(name, "DOT") == 0)
    {
        lineStyleFound = true;
    }
    else if (strcmp(name, "DOTDASH") == 0)
    {
        lineStyleFound = true;
    }
    if (lineStyleFound)
    {
        fprintf(stdout, "lineStyle is no longer a valid Histogram "
                       "attribute.\nIt's value is being ignored, please remove "
                       "it from your script.\n");
        return PyInt_FromLong(0L);
    }
    return Py_FindMethod(PyHistogramAttributes_methods, self, name);
}

int
PyHistogramAttributes_setattr(PyObject *self, char *name, PyObject *args)
{
    // Create a tuple to contain the arguments since all of the Set
    // functions expect a tuple.
    PyObject *tuple = PyTuple_New(1);
    PyTuple_SET_ITEM(tuple, 0, args);
    Py_INCREF(args);
    PyObject *obj = NULL;

    if(strcmp(name, "basedOn") == 0)
        obj = HistogramAttributes_SetBasedOn(self, tuple);
    else if(strcmp(name, "histogramType") == 0)
        obj = HistogramAttributes_SetHistogramType(self, tuple);
    else if(strcmp(name, "weightVariable") == 0)
        obj = HistogramAttributes_SetWeightVariable(self, tuple);
    else if(strcmp(name, "limitsMode") == 0)
        obj = HistogramAttributes_SetLimitsMode(self, tuple);
    else if(strcmp(name, "minFlag") == 0)
        obj = HistogramAttributes_SetMinFlag(self, tuple);
    else if(strcmp(name, "maxFlag") == 0)
        obj = HistogramAttributes_SetMaxFlag(self, tuple);
    else if(strcmp(name, "min") == 0)
        obj = HistogramAttributes_SetMin(self, tuple);
    else if(strcmp(name, "max") == 0)
        obj = HistogramAttributes_SetMax(self, tuple);
    else if(strcmp(name, "numBins") == 0)
        obj = HistogramAttributes_SetNumBins(self, tuple);
    else if(strcmp(name, "domain") == 0)
        obj = HistogramAttributes_SetDomain(self, tuple);
    else if(strcmp(name, "zone") == 0)
        obj = HistogramAttributes_SetZone(self, tuple);
    else if(strcmp(name, "useBinWidths") == 0)
        obj = HistogramAttributes_SetUseBinWidths(self, tuple);
    else if(strcmp(name, "outputType") == 0)
        obj = HistogramAttributes_SetOutputType(self, tuple);
    else if(strcmp(name, "lineWidth") == 0)
        obj = HistogramAttributes_SetLineWidth(self, tuple);
    else if(strcmp(name, "color") == 0)
        obj = HistogramAttributes_SetColor(self, tuple);
    else if(strcmp(name, "dataScale") == 0)
        obj = HistogramAttributes_SetDataScale(self, tuple);
    else if(strcmp(name, "binScale") == 0)
        obj = HistogramAttributes_SetBinScale(self, tuple);
    else if(strcmp(name, "normalizeHistogram") == 0)
        obj = HistogramAttributes_SetNormalizeHistogram(self, tuple);
    else if(strcmp(name, "computeAsCDF") == 0)
        obj = HistogramAttributes_SetComputeAsCDF(self, tuple);

    // Try and handle legacy fields
    if(obj == NULL)
    {
        if(strcmp(name, "lineStyle") == 0)
        {
            Py_INCREF(Py_None);
            obj = Py_None;
        }
    }
    if(obj != NULL)
        Py_DECREF(obj);

    Py_DECREF(tuple);
    if( obj == NULL)
        PyErr_Format(PyExc_RuntimeError, "Unable to set unknown attribute: '%s'", name);
    return (obj != NULL) ? 0 : -1;
}

static int
HistogramAttributes_print(PyObject *v, FILE *fp, int flags)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)v;
    fprintf(fp, "%s", PyHistogramAttributes_ToString(obj->data, "").c_str());
    return 0;
}

PyObject *
HistogramAttributes_str(PyObject *v)
{
    HistogramAttributesObject *obj = (HistogramAttributesObject *)v;
    return PyString_FromString(PyHistogramAttributes_ToString(obj->data,"").c_str());
}

//
// The doc string for the class.
//
#if PY_MAJOR_VERSION > 2 || (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION >= 5)
static const char *HistogramAttributes_Purpose = "Attributes for Histogram Plot";
#else
static char *HistogramAttributes_Purpose = "Attributes for Histogram Plot";
#endif

//
// The type description structure
//
static PyTypeObject HistogramAttributesType =
{
    //
    // Type header
    //
    PyObject_HEAD_INIT(&PyType_Type)
    0,                                   // ob_size
    "HistogramAttributes",                    // tp_name
    sizeof(HistogramAttributesObject),        // tp_basicsize
    0,                                   // tp_itemsize
    //
    // Standard methods
    //
    (destructor)HistogramAttributes_dealloc,  // tp_dealloc
    (printfunc)HistogramAttributes_print,     // tp_print
    (getattrfunc)PyHistogramAttributes_getattr, // tp_getattr
    (setattrfunc)PyHistogramAttributes_setattr, // tp_setattr
    (cmpfunc)HistogramAttributes_compare,     // tp_compare
    (reprfunc)0,                         // tp_repr
    //
    // Type categories
    //
    0,                                   // tp_as_number
    0,                                   // tp_as_sequence
    0,                                   // tp_as_mapping
    //
    // More methods
    //
    0,                                   // tp_hash
    0,                                   // tp_call
    (reprfunc)HistogramAttributes_str,        // tp_str
    0,                                   // tp_getattro
    0,                                   // tp_setattro
    0,                                   // tp_as_buffer
    Py_TPFLAGS_CHECKTYPES,               // tp_flags
    HistogramAttributes_Purpose,              // tp_doc
    0,                                   // tp_traverse
    0,                                   // tp_clear
    0,                                   // tp_richcompare
    0                                    // tp_weaklistoffset
};

//
// Helper functions for object allocation.
//

static HistogramAttributes *defaultAtts = 0;
static HistogramAttributes *currentAtts = 0;

static PyObject *
NewHistogramAttributes(int useCurrent)
{
    HistogramAttributesObject *newObject;
    newObject = PyObject_NEW(HistogramAttributesObject, &HistogramAttributesType);
    if(newObject == NULL)
        return NULL;
    if(useCurrent && currentAtts != 0)
        newObject->data = new HistogramAttributes(*currentAtts);
    else if(defaultAtts != 0)
        newObject->data = new HistogramAttributes(*defaultAtts);
    else
        newObject->data = new HistogramAttributes;
    newObject->owns = true;
    newObject->parent = 0;
    return (PyObject *)newObject;
}

static PyObject *
WrapHistogramAttributes(const HistogramAttributes *attr)
{
    HistogramAttributesObject *newObject;
    newObject = PyObject_NEW(HistogramAttributesObject, &HistogramAttributesType);
    if(newObject == NULL)
        return NULL;
    newObject->data = (HistogramAttributes *)attr;
    newObject->owns = false;
    newObject->parent = 0;
    return (PyObject *)newObject;
}

///////////////////////////////////////////////////////////////////////////////
//
// Interface that is exposed to the VisIt module.
//
///////////////////////////////////////////////////////////////////////////////

PyObject *
HistogramAttributes_new(PyObject *self, PyObject *args)
{
    int useCurrent = 0;
    if (!PyArg_ParseTuple(args, "i", &useCurrent))
    {
        if (!PyArg_ParseTuple(args, ""))
            return NULL;
        else
            PyErr_Clear();
    }

    return (PyObject *)NewHistogramAttributes(useCurrent);
}

//
// Plugin method table. These methods are added to the visitmodule's methods.
//
static PyMethodDef HistogramAttributesMethods[] = {
    {"HistogramAttributes", HistogramAttributes_new, METH_VARARGS},
    {NULL,      NULL}        /* Sentinel */
};

static Observer *HistogramAttributesObserver = 0;

std::string
PyHistogramAttributes_GetLogString()
{
    std::string s("HistogramAtts = HistogramAttributes()\n");
    if(currentAtts != 0)
        s += PyHistogramAttributes_ToString(currentAtts, "HistogramAtts.");
    return s;
}

static void
PyHistogramAttributes_CallLogRoutine(Subject *subj, void *data)
{
    typedef void (*logCallback)(const std::string &);
    logCallback cb = (logCallback)data;

    if(cb != 0)
    {
        std::string s("HistogramAtts = HistogramAttributes()\n");
        s += PyHistogramAttributes_ToString(currentAtts, "HistogramAtts.");
        cb(s);
    }
}

void
PyHistogramAttributes_StartUp(HistogramAttributes *subj, void *data)
{
    if(subj == 0)
        return;

    currentAtts = subj;
    PyHistogramAttributes_SetDefaults(subj);

    //
    // Create the observer that will be notified when the attributes change.
    //
    if(HistogramAttributesObserver == 0)
    {
        HistogramAttributesObserver = new ObserverToCallback(subj,
            PyHistogramAttributes_CallLogRoutine, (void *)data);
    }

}

void
PyHistogramAttributes_CloseDown()
{
    delete defaultAtts;
    defaultAtts = 0;
    delete HistogramAttributesObserver;
    HistogramAttributesObserver = 0;
}

PyMethodDef *
PyHistogramAttributes_GetMethodTable(int *nMethods)
{
    *nMethods = 1;
    return HistogramAttributesMethods;
}

bool
PyHistogramAttributes_Check(PyObject *obj)
{
    return (obj->ob_type == &HistogramAttributesType);
}

HistogramAttributes *
PyHistogramAttributes_FromPyObject(PyObject *obj)
{
    HistogramAttributesObject *obj2 = (HistogramAttributesObject *)obj;
    return obj2->data;
}

PyObject *
PyHistogramAttributes_New()
{
    return NewHistogramAttributes(0);
}

PyObject *
PyHistogramAttributes_Wrap(const HistogramAttributes *attr)
{
    return WrapHistogramAttributes(attr);
}

void
PyHistogramAttributes_SetParent(PyObject *obj, PyObject *parent)
{
    HistogramAttributesObject *obj2 = (HistogramAttributesObject *)obj;
    obj2->parent = parent;
}

void
PyHistogramAttributes_SetDefaults(const HistogramAttributes *atts)
{
    if(defaultAtts)
        delete defaultAtts;

    defaultAtts = new HistogramAttributes(*atts);
}

