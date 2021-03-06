// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef RESAMPLEATTRIBUTES_H
#define RESAMPLEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: ResampleAttributes
//
// Purpose:
//    Atts for Resample operator
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class ResampleAttributes : public AttributeSubject
{
public:
    enum TieResolver
    {
        random,
        largest,
        smallest
    };

    // These constructors are for objects of this class
    ResampleAttributes();
    ResampleAttributes(const ResampleAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    ResampleAttributes(private_tmfs_t tmfs);
    ResampleAttributes(const ResampleAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~ResampleAttributes();

    virtual ResampleAttributes& operator = (const ResampleAttributes &obj);
    virtual bool operator == (const ResampleAttributes &obj) const;
    virtual bool operator != (const ResampleAttributes &obj) const;
private:
    void Init();
    void Copy(const ResampleAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectTieResolverVariable();

    // Property setting methods
    void SetUseExtents(bool useExtents_);
    void SetStartX(double startX_);
    void SetEndX(double endX_);
    void SetSamplesX(int samplesX_);
    void SetStartY(double startY_);
    void SetEndY(double endY_);
    void SetSamplesY(int samplesY_);
    void SetIs3D(bool is3D_);
    void SetStartZ(double startZ_);
    void SetEndZ(double endZ_);
    void SetSamplesZ(int samplesZ_);
    void SetTieResolver(TieResolver tieResolver_);
    void SetTieResolverVariable(const std::string &tieResolverVariable_);
    void SetDefaultValue(double defaultValue_);
    void SetDistributedResample(bool distributedResample_);
    void SetCellCenteredOutput(bool cellCenteredOutput_);

    // Property getting methods
    bool              GetUseExtents() const;
    double            GetStartX() const;
    double            GetEndX() const;
    int               GetSamplesX() const;
    double            GetStartY() const;
    double            GetEndY() const;
    int               GetSamplesY() const;
    bool              GetIs3D() const;
    double            GetStartZ() const;
    double            GetEndZ() const;
    int               GetSamplesZ() const;
    TieResolver       GetTieResolver() const;
    const std::string &GetTieResolverVariable() const;
          std::string &GetTieResolverVariable();
    double            GetDefaultValue() const;
    bool              GetDistributedResample() const;
    bool              GetCellCenteredOutput() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string TieResolver_ToString(TieResolver);
    static bool TieResolver_FromString(const std::string &, TieResolver &);
protected:
    static std::string TieResolver_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    virtual bool SetValue(const std::string &name, const double &value);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_useExtents = 0,
        ID_startX,
        ID_endX,
        ID_samplesX,
        ID_startY,
        ID_endY,
        ID_samplesY,
        ID_is3D,
        ID_startZ,
        ID_endZ,
        ID_samplesZ,
        ID_tieResolver,
        ID_tieResolverVariable,
        ID_defaultValue,
        ID_distributedResample,
        ID_cellCenteredOutput,
        ID__LAST
    };

private:
    bool        useExtents;
    double      startX;
    double      endX;
    int         samplesX;
    double      startY;
    double      endY;
    int         samplesY;
    bool        is3D;
    double      startZ;
    double      endZ;
    int         samplesZ;
    int         tieResolver;
    std::string tieResolverVariable;
    double      defaultValue;
    bool        distributedResample;
    bool        cellCenteredOutput;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define RESAMPLEATTRIBUTES_TMFS "bddiddibddiisdbb"

#endif
