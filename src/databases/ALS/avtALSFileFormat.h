// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                            avtALSFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_ALS_FILE_FORMAT_H
#define AVT_ALS_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>

#include <vector>
#include <hdf5.h>

// ****************************************************************************
//  Class: avtALSFileFormat
//
//  Purpose:
//      Reads in ALS files as a plugin to VisIt.
//
//  Programmer: hari -- generated by xml2avt
//  Creation:   Mon Apr 14 16:36:42 PST 2014
//
// ****************************************************************************

class avtALSFileFormat : public avtSTMDFileFormat
{
  public:
                       avtALSFileFormat(const char *);
    virtual           ~avtALSFileFormat() {;}

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
    // virtual void      *GetAuxiliaryData(const char *var, int domain,
    //                                     const char *type, void *args, 
    //                                     DestructorFunction &);
    //

    //
    // If you know the cycle number, overload this function.
    // Otherwise, VisIt will make up a reasonable one for you.
    //
    // virtual int         GetCycle(void);
    //

    virtual const char    *GetType(void)   { return "ALS"; }
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

    void SetGroupName(const std::string& name) { m_groupName = name; }

  protected:
    enum ALSDataType {STANDARD, TOMO, OTHER};

    ALSDataType m_dataType;
    // DATA MEMBERS
    std::string            m_filename;
    size_t                 m_width, m_height, m_slices;
    std::string            m_groupName;
    bool                   m_initialized;
    std::vector<std::string> m_datasetNames;
    int                    m_ver;
    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);

  private:
    void Initialize();
    bool InitializeHeader();
    bool InitializeStandardHeader(hid_t file);
    bool InitializeTomoHeader(hid_t file);
    vtkDataArray* GetDataSet();
    vtkDataArray* GetTomoDataSet(hid_t file);
    vtkDataArray* GetStandardDataSet(hid_t file);
};


#endif
