// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                            avtMFEMFileFormat.C                           //
// ************************************************************************* //

#include <avtMFEMFileFormat.h>

#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>

#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>

#include <avtDatabaseMetaData.h>

#include <DBOptionsAttributes.h>
#include <DebugStream.h>
#include <Expression.h>
#include <FileFunctions.h>

#include <InvalidFilesException.h>
#include <InvalidVariableException.h>

#include <avtResolutionSelection.h>

#include <StringHelpers.h>
#include <visit_gzstream.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkHexahedron.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkPoints.h>
#include "mfem.hpp"

#include <JSONRoot.h>

#include <zlib.h>

#ifdef _WIN32
#define strncasecmp _strnicmp
static unsigned long long
v_strtoull(const char *__restrict str, char **__restrict endptr, int base)
{
    return (unsigned long long)strtoul(str, endptr, base);
}
#else
#define v_strtoull strtoull
#endif

using     std::string;
using     std::ostringstream;
using     std::vector;

using namespace mfem;

// ****************************************************************************
//  Method: avtMFEMFileFormat constructor
//
//  Programmer: harrison37 -- generated by xml2avt
//  Creation:   Fri May 23 15:16:20 PST 2014
//
// ****************************************************************************

avtMFEMFileFormat::avtMFEMFileFormat(const char *filename)
    : avtSTMDFileFormat(&filename, 1)
{
    selectedLOD = 0;
    root        = NULL;
}

// ****************************************************************************
//  Method: avtMFEMFileFormat destructor
//
//  Programmer: harrison37 -- generated by xml2avt
//  Creation:   Fri May 23 15:16:20 PST 2014
//
// ****************************************************************************

avtMFEMFileFormat::~avtMFEMFileFormat()
{
    FreeUpResources();
}


// ****************************************************************************
//  Method: avtMFEMFileFormat::FreeUpResources
//
//  Purpose:
//      When VisIt is done focusing on a particular timestep, it asks that
//      timestep to free up any resources (memory, file descriptors) that
//      it has associated with it.  This method is the mechanism for doing
//      that.
//
//  Programmer: harrison37 -- generated by xml2avt
//  Creation:   Fri May 23 15:16:20 PST 2014
//
// ****************************************************************************

void
avtMFEMFileFormat::ActivateTimestep(void)
{
    avtSTMDFileFormat::ActivateTimestep();
    if(root)
        FreeUpResources();
    std::string root_file(GetFilename());
    root = new JSONRoot(root_file);
}


void
avtMFEMFileFormat::FreeUpResources(void)
{
    if(root)
    {
        delete root;
        root = NULL;
    }
}


// ****************************************************************************
//  Method: avtMFEMFileFormat::BuildCatFileMap
//
//  Purpose: Read .mfem_cat file header for offsets/sizes of sub-files
//
//  Mark C. Miller, Mon Dec 11 15:48:47 PST 2017
//
//  Mark C. Miller, Wed Dec 13 16:37:35 PST 2017
//  Add support for compression
// ****************************************************************************
void
avtMFEMFileFormat::BuildCatFileMap(string const &cat_path)
{
    if (catFileMap.size())
        return;

    ifstream catfile(cat_path.c_str());
    if (!catfile)
    {
        debug5 << "Failed to open mfem_cat file, \"" << cat_path << "\"";
        return;
    }

    string line;
    std::getline(catfile, line);
    size_t hdrsz = (size_t) v_strtoull(&line[0], 0, 10);
    size_t zip = strchr(&line[0], 'z')==0?(size_t)0:(size_t)1;
    catFileMap["@header_size@"] = std::pair<size_t,size_t>(hdrsz,hdrsz);
    catFileMap["@compressed@"] = std::pair<size_t,size_t>(zip,zip);
    debug5 << "Processing mfem_cat file header..." << endl;
    debug5 << "    header size = " << hdrsz << endl;
    debug5 << "    compressed  = " << zip << endl;
    while (catfile.tellg() < hdrsz)
    {
        std::getline(catfile, line);
        size_t offat = line.find_last_of(' ');
        size_t sizat = line.find_last_of(' ', offat-1);
        size_t off = (size_t) v_strtoull(&line[offat+1], 0, 10);
        size_t siz = (size_t) v_strtoull(&line[sizat+1], 0, 10);
        line.resize(sizat);
        debug5 << "    key=\"" << line << "\", size=" << siz << ", off=" << off << endl;
        catFileMap[line] = std::pair<size_t,size_t>(siz,off);
    }
}

// ****************************************************************************
//  Method: avtMFEMFileFormat::FetchDataFromCatFile
//
//  Purpose: Read mfmem data from an mfem_cat file
//
//  Mark C. Miller, Mon Dec 11 15:48:47 PST 2017
// ****************************************************************************
void
avtMFEMFileFormat::FetchDataFromCatFile(string const &cat_path, string const &obj_path,
    std::istringstream &istr)
{

    BuildCatFileMap(cat_path);

    string line;

    // Indicate error if we are unable to locate the sub-file
    string obj_base = FileFunctions::Basename(obj_path);
    if (catFileMap.find(obj_base) == catFileMap.end())
        return istr.setstate(std::ios::failbit);

    // Look up the sub-file size and offset and header_size
    size_t size = catFileMap[obj_base].first;
    size_t offset = catFileMap[obj_base].second;
    size_t header_size = catFileMap["@header_size@"].first; // special key for hdrsz
    bool zip = catFileMap["@compressed@"].first != 0; // special key for hdrsz

    debug5 << "Fetching " << (zip?"compressed":"") << " data of size " 
           << size << " at offset " << offset+header_size << endl;

    ifstream catfile(cat_path.c_str());
    if (!catfile)
        return istr.setstate(std::ios::failbit);

    string data;
    data.resize(size+1);
    catfile.seekg(offset+header_size);
    catfile.read(&data[0], size);

    if (!zip)
    {
        istr.str(data);
        return;
    }

    z_stream zstrm;

    zstrm.zalloc = Z_NULL;
    zstrm.zfree = Z_NULL;
    zstrm.opaque = Z_NULL;
    zstrm.next_in = Z_NULL;
    zstrm.avail_in = 0;

    if (inflateInit2(&zstrm, 15+16) != Z_OK)
        return istr.setstate(std::ios::failbit);

    string newdata;
    zstrm.avail_in = size;
    zstrm.next_in = (Bytef*) &data[0];
    int sum = 0;
    while (true)
    {
        int const bufsiz = 10*1024;
        char buf[bufsiz];

        zstrm.avail_out = bufsiz;
        zstrm.next_out = (Bytef*)buf;
        int zstatus = inflate(&zstrm, Z_NO_FLUSH);
        if (zstatus == Z_STREAM_ERROR || zstatus == Z_NEED_DICT ||
            zstatus == Z_DATA_ERROR || zstatus == Z_MEM_ERROR)
        {
            inflateEnd(&zstrm);
            return istr.setstate(std::ios::failbit);
        }
        int have = bufsiz - zstrm.avail_out;
        newdata.append(buf, have);
        sum += have;
        if (zstatus == Z_STREAM_END)
            break;
    }

    if (inflateEnd(&zstrm) != Z_OK)
        return istr.setstate(std::ios::failbit);

    debug5 << "Decompressed " << size << " bytes into " << sum << " bytes." << endl;
    istr.str(newdata);
}

// ****************************************************************************
//  Method: avtMFEMFileFormat::PopulateDatabaseMetaData
//
//  Purpose:
//      This database meta-data object is like a table of contents for the
//      file.  By populating it, you are telling the rest of VisIt what
//      information it can request from you.
//
//  Programmer: harrison37 -- generated by xml2avt
//  Creation:   Fri May 23 15:16:20 PST 2014
//
//
//  Modifications:
//   Cyrus Harrison, Wed Sep 24 10:47:00 PDT 2014
//   Move abs path logic into JSONRoot.
//
//   Mark C. Miller, Tue Sep 20 18:12:29 PDT 2016
//   Add expressions
// ****************************************************************************
void
avtMFEMFileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md)
{    
    ///
    /// Open the root file
    ///
    std::string root_file(GetFilename());

    JSONRoot root_md(root_file);
    vector<string>dset_names;
    root_md.DataSets(dset_names);

    // enumerate datasets)
    selectedLOD = 0;
    for(int i=0;i<(int)dset_names.size();i++)
    {
        JSONRootDataSet &dset =  root_md.DataSet(dset_names[i]);
        int nblocks      = dset.NumberOfDomains();
        int spatial_dim  = atoi(dset.Mesh().Tag("spatial_dim").c_str());
        int topo_dim     = atoi(dset.Mesh().Tag("topo_dim").c_str());
        int block_origin = 0;
        double *extents = NULL;

        //
        // Add the meta-data object for the current mesh:
        //
        AddMeshToMetaData(md,
                          dset_names[i].c_str(),
                          AVT_UNSTRUCTURED_MESH,
                          extents,
                          nblocks,
                          block_origin,
                          spatial_dim, topo_dim);

        md->GetMeshes(i).LODs = atoi(dset.Mesh().Tag("max_lods").c_str());

        // Add builtin mfem fields related to the mesh:

        AddScalarVarToMetaData(md,
                               "element_coloring",
                               dset_names[i].c_str(),
                               AVT_ZONECENT);
        AddScalarVarToMetaData(md,
                               "element_attribute",
                               dset_names[i].c_str(),
                               AVT_ZONECENT);
    
    
        /// add mesh variables
        vector<string>field_names;
        dset.Fields(field_names);
        for(size_t j=0;j<field_names.size();j++)
        {
            JSONRootEntry &field = dset.Field(field_names[j]);
            std::string slod = field.Tag("lod");
            int ilod = std::min(md->GetMeshes(i).LODs,atoi(slod.c_str()));
            selectedLOD = std::max(selectedLOD,ilod);
            std::string f_assoc = field.Tag("assoc");

            if(f_assoc == "elements")
            {
                if(!field.HasTag("comps") || field.Tag("comps") == "1")
                {
                    AddScalarVarToMetaData(md,
                                           field_names[j].c_str(),
                                           dset_names[i].c_str(),
                                           AVT_ZONECENT);
                }
                else if(field.Tag("comps") == "2")
                {
                    AddVectorVarToMetaData(md,
                                          field_names[j].c_str(),
                                          dset_names[i].c_str(),
                                          AVT_ZONECENT,2);
                }
                else if(field.Tag("comps") == "3")
                {
                    AddVectorVarToMetaData(md,
                                          field_names[j].c_str(),
                                          dset_names[i].c_str(),
                                          AVT_ZONECENT,3);
                }


            }
            else if(f_assoc == "nodes")
            {
                if(field.Tag("comps") == "1")
                {
                    AddScalarVarToMetaData(md,
                                           field_names[j].c_str(),
                                           dset_names[i].c_str(),
                                           AVT_NODECENT);
                }
                else if(field.Tag("comps") == "2")
                {
                    AddVectorVarToMetaData(md,
                                          field_names[j].c_str(),
                                          dset_names[i].c_str(),
                                          AVT_NODECENT,2);
                }
                else if(field.Tag("comps") == "3")
                {
                    AddVectorVarToMetaData(md,
                                          field_names[j].c_str(),
                                          dset_names[i].c_str(),
                                          AVT_NODECENT,3);
                }
            }
        }
    }

    vector<string>expr_names;
    root_md.Expressions(expr_names);
    for(size_t i = 0; i < expr_names.size(); i++)
    {
        JSONRootExpr &jexpr = root_md.Expression(expr_names[i]);
        
        Expression::ExprType vartype = Expression::Unknown;
        if      (!strncasecmp(jexpr.Type().c_str(), "scalar", 6))
            vartype = Expression::ScalarMeshVar;
        else if (!strncasecmp(jexpr.Type().c_str(), "vector", 6))
            vartype = Expression::VectorMeshVar;
        else if (!strncasecmp(jexpr.Type().c_str(), "tensor", 6))
            vartype = Expression::TensorMeshVar;
        else if (!strncasecmp(jexpr.Type().c_str(), "array", 5))
            vartype = Expression::ArrayMeshVar;
        else if (!strncasecmp(jexpr.Type().c_str(), "material", 8))
            vartype = Expression::Material;
        else if (!strncasecmp(jexpr.Type().c_str(), "species", 7))
            vartype = Expression::Species;
        else
        {
            debug5 << "Warning: unknown expression type \"" << jexpr.Type() << "\" for expression "
                   << "\"" << expr_names[i] << "\"...skipping it." << endl;
            continue;
        }

        Expression expr;
        expr.SetName(expr_names[i]);
        expr.SetDefinition(jexpr.Defn());
        expr.SetType(vartype);
        md->AddExpression(&expr);
    }
}

// ****************************************************************************
//  Method: avtMFEMFileFormat::GetCycle
//
//  Purpose:
//      Returns if we have the current cycle.
//
//  Programmer: Cyrus Harrison
//  Creation:   Wed Oct 15 10:52:22 PDT 2014
//
// ****************************************************************************
int
avtMFEMFileFormat::GetCycle() 
{
    // VisIt doesn't support diff times / cycles for meshes
    // we loop over all meshes to see if any have valid cycle info
    if(root)
    {
        std::vector<std::string> dset_names;
        root->DataSets(dset_names);
        // enumerate datasets
        bool not_found = true;
        for(size_t i=0; i<dset_names.size() && not_found ;i++)
        {
            JSONRootDataSet &dset =  root->DataSet(dset_names[i]);
            if(dset.HasCycle())
            {
                return dset.Cycle();
            }
        }
    }
    
    return avtFileFormat::INVALID_CYCLE;
}


// ****************************************************************************
//  Method: avtMFEMFileFormat::GetTime
//
//  Purpose:
//      Returns if we have the current cycle.
//
//  Programmer: Cyrus Harrison
//  Creation:   Wed Oct 15 10:52:22 PDT 2014
//
// ****************************************************************************
double
avtMFEMFileFormat::GetTime() 
{
    // VisIt doesn't support diff times / cycles for meshes
    // we loop over all meshes to see if any have valid time info
    if(root)
    {
        std::vector<std::string> dset_names;
        root->DataSets(dset_names);
        // enumerate datasets
        bool not_found = true;
        for(size_t i=0; i<dset_names.size() && not_found ;i++)
        {
            JSONRootDataSet &dset =  root->DataSet(dset_names[i]);
            if(dset.HasTime())
            {
                return dset.Time();
            }
        }
    }
    
    return avtFileFormat::INVALID_TIME;
}


// ****************************************************************************
//  Method: avtMFEMFileFormat::GetMesh
//
//  Purpose:
//      Gets the mesh associated with this file.  The mesh is returned as a
//      derived type of vtkDataSet (ie vtkRectilinearGrid, vtkStructuredGrid,
//      vtkUnstructuredGrid, etc).
//
//  Arguments:
//      domain      The index of the domain.  If there are NDomains, this
//                  value is guaranteed to be between 0 and NDomains-1,
//                  regardless of block origin.
//      meshname    The name of the mesh of interest.  This can be ignored if
//                  there is only one mesh.
//
//  Programmer: harrison37 -- generated by xml2avt
//  Creation:   Fri May 23 15:16:20 PST 2014
//
// ****************************************************************************
vtkDataSet *
avtMFEMFileFormat::GetMesh(int domain, const char *meshname)
{
    return GetRefinedMesh(string(meshname),domain,selectedLOD+1);
}



// ****************************************************************************
//  Method: avtMFEMFileFormat::GetVar
//
//  Purpose:
//      Gets a scalar variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      domain     The index of the domain.  If there are NDomains, this
//                 value is guaranteed to be between 0 and NDomains-1,
//                 regardless of block origin.
//      varname    The name of the variable requested.
//
//  Programmer: harrison37 -- generated by xml2avt
//  Creation:   Fri May 23 15:16:20 PST 2014
//
// ****************************************************************************
vtkDataArray *
avtMFEMFileFormat::GetVar(int domain, const char *varname)
{   
   return GetRefinedVar(string(varname),domain,selectedLOD+1);
}

// ****************************************************************************
//  Method: avtMFEMFileFormat::GetVectorVar
//
//  Purpose:
//      Gets a vector variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      domain     The index of the domain.  If there are NDomains, this
//                 value is guaranteed to be between 0 and NDomains-1,
//                 regardless of block origin.
//      varname    The name of the variable requested.
//
//  Programmer: harrison37 -- generated by xml2avt
//  Creation:   Fri May 23 15:16:20 PST 2014
//
// ****************************************************************************

vtkDataArray *
avtMFEMFileFormat::GetVectorVar(int domain, const char *varname)
{
    return GetRefinedVar(string(varname),domain,selectedLOD+1);
}


// ****************************************************************************
//  Method: avtMFEMFileFormat::FetchMesh
//
//  Purpose: Returns a new instance of the mfem mesh, given a domain id and 
//           mesh name
//
//  Arguments:
//    mesh_name: string with desired mesh name
//    domain:    domain id
//
//
//  Programmer: Cyrus Harrison
//  Creation:   Sat Jul  5 11:38:31 PDT 2014
//
// Modifications:
//   Cyrus Harrison, Wed Jun  1 08:46:12 PDT 2016
//   Change MFEM Mesh constructor call to resolve coordinate system issue
//   (See: http://visitbugs.ornl.gov/issues/2578)
//
//   Cyrus Harrison, Mon Aug 22 20:00:57 PDT 2016
//   Additional change to MFEM Mesh constructor call to resolve 
//   coordinate system issue
//
//   Mark C. Miller, Mon Dec 11 15:49:34 PST 2017
//   Add support for mfem_cat file
// ****************************************************************************
Mesh *
avtMFEMFileFormat::FetchMesh(const std::string &mesh_name,int domain)
{
    if(root == NULL)
    {
        //failed to open mesh file
        ostringstream msg;
        msg << "Failed to open MFEM mesh:"
            << " root metadata is missing!";
        EXCEPTION1(InvalidFilesException, msg.str());
    }

    string mesh_path = root->DataSet(mesh_name).Mesh().Path().Expand(domain);
    string cat_path = root->DataSet(mesh_name).CatPath().Get();

    if (cat_path != "")
    {
        std::istringstream imeshstr;
        FetchDataFromCatFile(cat_path, mesh_path, imeshstr); 

        // failed to get to mesh data from cat file
        if (imeshstr)
            return new Mesh(imeshstr, 1, 0, false);
    }

    visit_ifstream imesh(mesh_path.c_str());
    if(imesh().fail())
    {
        //failed to open mesh file
        ostringstream msg;
        msg << "Failed to open MFEM mesh:"
            << " mesh name: " << mesh_name 
            << " domain: "    << domain
            << " mesh path: \"" << mesh_path << "\"";
        if (cat_path != "")
            msg << " cat path: \"" << cat_path << "\"";
        EXCEPTION1(InvalidFilesException, msg.str());
    }
   
    return new Mesh(imesh(), 1, 0, false);
}

// ****************************************************************************
//  Method: avtMFEMFileFormat::GetRefinedMesh
//
//  Purpose: 
//    Constructs a vtkUnstructuredGrid that contains a refined mfem mesh.
//
//  Arguments:
//    mesh_name: string with desired mesh name
//    domain:    domain id
//    lod:       number of refinement steps 
//
//  Programmer: Cyrus Harrison
//  Creation:   Sat Jul  5 11:38:31 PDT 2014
//
// ****************************************************************************
vtkDataSet *
avtMFEMFileFormat::GetRefinedMesh(const std::string &mesh_name, int domain, int lod)
{
     // get base mesh
    Mesh *mesh = FetchMesh(mesh_name,domain);
     
    // create output objects
    vtkUnstructuredGrid *res_ds  = vtkUnstructuredGrid::New(); 
    vtkPoints           *res_pts = vtkPoints::New();
   
    int npts=0;
    int neles=0;
    
    RefinedGeometry *refined_geo;
    DenseMatrix      pmat;
   
    //
    // find the # of output points and cells at the selected level of 
    // refinement (we may want to cache this...)
    //
    for (int i = 0; i < mesh->GetNE(); i++)
    {
        int geom = mesh->GetElementBaseGeometry(i);
        int ele_nverts = Geometries.GetVertices(geom)->GetNPoints();
        refined_geo = GlobGeometryRefiner.Refine(geom, lod, 1);
        npts  += refined_geo->RefPts.GetNPoints();
        neles += refined_geo->RefGeoms.Size() / ele_nverts;
    }

    // create the points for the refined topoloy   
    res_pts->Allocate(npts);
    res_pts->SetNumberOfPoints((vtkIdType) npts);
   
    // create the points for the refined topoloy
    int pt_idx=0;
    for (int i = 0; i < mesh->GetNE(); i++)
    {
        int geom = mesh->GetElementBaseGeometry(i);
        refined_geo = GlobGeometryRefiner.Refine(geom, lod, 1);
        // refined points
        mesh->GetElementTransformation(i)->Transform(refined_geo->RefPts, pmat);
        for (int j = 0; j < pmat.Width(); j++)
        {
            double pt_vals[3]={0.0,0.0,0.0};
            pt_vals[0] = pmat(0,j);
            if (pmat.Height() > 1)
                pt_vals[1] = pmat(1,j);
            if (pmat.Height() > 2)
                pt_vals[2] = pmat(2,j);
            res_pts->InsertPoint(pt_idx,pt_vals);
            pt_idx++;
        }
    }
    
    res_ds->SetPoints(res_pts);
    res_pts->Delete();  
    // create the cells for the refined topology   
    res_ds->Allocate(neles);
    
    pt_idx=0;
    for (int i = 0; i <  mesh->GetNE(); i++)
    {
        int geom       = mesh->GetElementBaseGeometry(i);
        int ele_nverts = Geometries.GetVertices(geom)->GetNPoints();
        refined_geo    = GlobGeometryRefiner.Refine(geom, lod, 1);

        Array<int> &rg_idxs = refined_geo->RefGeoms;

        vtkCell *ele_cell = NULL;
        // rg_idxs contains all of the verts for the refined elements
        for (int j = 0; j < rg_idxs.Size(); )
        {
            switch (geom)
            {
                case Geometry::SEGMENT:      ele_cell = vtkLine::New();       break;
                case Geometry::TRIANGLE:     ele_cell = vtkTriangle::New();   break;
                case Geometry::SQUARE:       ele_cell = vtkQuad::New();       break;
                case Geometry::TETRAHEDRON:  ele_cell = vtkTetra::New();      break;
                case Geometry::CUBE:         ele_cell = vtkHexahedron::New(); break;
            }
            // the are ele_nverts for each refined element
            for (int k = 0; k < ele_nverts; k++, j++)
                ele_cell->GetPointIds()->SetId(k,pt_idx + rg_idxs[j]);

            res_ds->InsertNextCell(ele_cell->GetCellType(),
                                   ele_cell->GetPointIds());
            ele_cell->Delete();
        }

        pt_idx += refined_geo->RefPts.GetNPoints();
   }
   
   delete mesh;
       
   return res_ds;
}

// ****************************************************************************
//  Method: avtMFEMFileFormat::GetRefinedVar
//
//  Purpose: 
//   Constructs a vtkDataArray that contains a refined mfem mesh field variable.
//
//  Arguments:
//   var_name:  string with desired variable name
//   domain:    domain id
//   lod:       number of refinement steps 
//
//  Programmer: Cyrus Harrison
//  Creation:   Sat Jul  5 11:38:31 PDT 2014
//
//  Mark C. Miller, Mon Dec 11 15:49:34 PST 2017
//  Add support for mfem_cat file
// ****************************************************************************
vtkDataArray *
avtMFEMFileFormat::GetRefinedVar(const std::string &var_name,
                                 int domain,
                                 int lod)
{
    if(root == NULL)
    {
        //failed to open mesh file
        ostringstream msg;
        msg << "Failed to fetch MFEM grid function:"
            << " root metadata is missing!";
        EXCEPTION1(InvalidFilesException, msg.str());
    }

    // find mesh for var
    string mesh_name="main";
    vector<string> dset_names;
    root->DataSets(dset_names);
    for(size_t i=0;i<dset_names.size();i++)
    {
        JSONRootDataSet &dset = root->DataSet(dset_names[i]);
        vector<string> field_names;
        dset.Fields(field_names);
        for(size_t j=0;j<field_names.size();j++)
        {
            if(field_names[j] == var_name)
                mesh_name = dset_names[i];
        }
    }

    // check for special vars
    if(var_name == "element_coloring")
        return GetRefinedElementColoring(mesh_name,domain,lod);
    else if(var_name == "element_attribute") // handle with materials in the future?
        return GetRefinedElementAttribute(mesh_name,domain,lod);

    // get base mesh
    Mesh *mesh = FetchMesh(mesh_name,domain);

    JSONRootEntry &field = root->DataSet(mesh_name).Field(var_name);
    string field_path = field.Path().Expand(domain);
    bool var_is_nodal = field.Tag("assoc") == "nodes";
    int  ncomps       = atoi(field.Tag("comps").c_str());
    string cat_path = root->DataSet(mesh_name).CatPath().Get();

    GridFunction *gf = 0;
    if (cat_path != "")
    {
        std::istringstream igfstr;
        FetchDataFromCatFile(cat_path, field_path, igfstr); 

        if (igfstr)
            gf = new GridFunction(mesh,igfstr);   
    }

    if (!gf)
    {
        visit_ifstream igf(field_path.c_str());
        if (igf().fail())
        {
            //failed to open gf file
            ostringstream msg;
            msg << "Failed to open MFEM grid function: "
                << " field name: \""       << mesh_name << "\""
                << " domain: "             << domain
                << " grid function path: \"" << field_path << "\"";
            if (cat_path != "")
                msg << " cat path: \"" << cat_path << "\"";

            EXCEPTION1(InvalidFilesException, msg.str());
        }
        gf = new GridFunction(mesh,igf());   
    }

    int npts=0;
    int neles=0;
    
    RefinedGeometry *refined_geo;
    Vector           scalar_vals;
    DenseMatrix      vec_vals;
    DenseMatrix      pmat;

    //
    // find the # of output points and cells at the selected level of 
    // refinement (we may want to cache this...)
    //
    for (int i = 0; i < mesh->GetNE(); i++)
    {
        int geom = mesh->GetElementBaseGeometry(i);
        int ele_nverts = Geometries.GetVertices(geom)->GetNPoints();
        refined_geo    = GlobGeometryRefiner.Refine(geom, lod, 1);
        npts  += refined_geo->RefPts.GetNPoints();
        neles += refined_geo->RefGeoms.Size() / ele_nverts;
    }
    
    vtkFloatArray *rv = vtkFloatArray::New();
    if(ncomps == 2)
        rv->SetNumberOfComponents(3);
    else
        rv->SetNumberOfComponents(ncomps);
    if(var_is_nodal)
        rv->SetNumberOfTuples(npts);
    else
        rv->SetNumberOfTuples(neles);

    double tuple_vals[9];
    int ref_idx=0;
    for (int i = 0; i <  mesh->GetNE(); i++)
    {
        int geom       = mesh->GetElementBaseGeometry(i);
        refined_geo    = GlobGeometryRefiner.Refine(geom, lod, 1);
        if(ncomps == 1)
        {
            gf->GetValues(i, refined_geo->RefPts, scalar_vals, pmat);
            for (int j = 0; j < scalar_vals.Size();j++)
            {
                tuple_vals[0] = scalar_vals(j);
                rv->SetTuple(ref_idx, tuple_vals); 
                ref_idx++;
            }
        }
        else
        {
            gf->GetVectorValues(i, refined_geo->RefPts, vec_vals, pmat);
            for (int j = 0; j < vec_vals.Width(); j++)
            {
                tuple_vals[2] = 0.0;
                tuple_vals[0] = vec_vals(0,j);
                tuple_vals[1] = vec_vals(1,j);
                if (vec_vals.Height() > 2)
                    tuple_vals[2] = vec_vals(2,j);
                rv->SetTuple(ref_idx, tuple_vals); 
                ref_idx++;
            }
        }
    }
    
    delete gf;
    delete mesh;

    return rv;
}

// ****************************************************************************
//  Method: avtMFEMFileFormat::GetRefinedElementColoring
//
//  Purpose:
//   Constructs a vtkDataArray that contains coloring that refects the orignal
//   finite elements of a mfem mesh.
//
//  Arguments:
//   var_name:  string with desired mesh name
//   domain:    domain id
//   lod:       number of refinement steps 
//
//  Programmer: Cyrus Harrison
//  Creation:   Sat Jul  5 11:38:31 PDT 2014
//
// Modifications:
//   Cyrus Harrison, Tue May 23 10:12:52 PDT 2017
//   Seed rng with domain id for predictable coloring results
//   (See: http://visitbugs.ornl.gov/issues/2747)
//
// ****************************************************************************
vtkDataArray *
avtMFEMFileFormat::GetRefinedElementColoring(const std::string &mesh_name,
                                             int domain,
                                             int lod)
{
    //
    // fetch the base mfem mesh
    //
    Mesh *mesh = FetchMesh(mesh_name,domain);
    int npts=0;
    int neles=0;
    
    RefinedGeometry *refined_geo;
    Array<int>       coloring;
    
    //
    // find the # of output points and cells at the selected level of 
    // refinement (we may want to cache this...)
    //
    for (int i = 0; i < mesh->GetNE(); i++)
    {
        int geom = mesh->GetElementBaseGeometry(i);
        int ele_nverts = Geometries.GetVertices(geom)->GetNPoints();
        refined_geo    = GlobGeometryRefiner.Refine(geom, lod, 1);
        npts  += refined_geo->RefPts.GetNPoints();
        neles += refined_geo->RefGeoms.Size() / ele_nverts;
    }

    vtkFloatArray *rv = vtkFloatArray::New();
    rv->SetNumberOfComponents(1);
    rv->SetNumberOfTuples(neles);

    //
    // Use mfem's mesh coloring alog
    //

    // seed using domain id for predictable results
    srand(domain);

    double a = double(rand()) / (double(RAND_MAX) + 1.);

    int el0 = (int)floor(a * mesh->GetNE());
    mesh->GetElementColoring(coloring, el0);
    int ref_idx=0;
    // set output array value from generated coloring
    for (int i = 0; i < mesh->GetNE(); i++)
    {
        int geom = mesh->GetElementBaseGeometry(i);
        int nv = Geometries.GetVertices(geom)->GetNPoints();
        refined_geo= GlobGeometryRefiner.Refine(geom, lod, 1);
        for (int j = 0; j < refined_geo->RefGeoms.Size(); j += nv)
        {
             rv->SetTuple1(ref_idx,coloring[i]+1);
             ref_idx++;
        }
   }
   
   delete mesh;
   
   return rv;
}


// ****************************************************************************
//  Method: avtMFEMFileFormat::GetRefinedElementAttribute
//
//  Purpose:
//   Constructs a vtkDataArray that contains the refined "attribute" value 
//   for finite elements in a mfem mesh.
//
//  Arguments:
//   var_name:  string with desired mesh name
//   domain:    domain id
//   lod:       number of refinement steps 
//
//  Programmer: Cyrus Harrison
//  Creation:   Sat Jul  5 11:38:31 PDT 2014
//
// ****************************************************************************
vtkDataArray *
avtMFEMFileFormat::GetRefinedElementAttribute(const std::string &mesh_name, 
                                              int domain, 
                                              int lod)
{
    //
    // fetch the base mfem mesh
    //
    Mesh *mesh = FetchMesh(mesh_name,domain);
    int npts=0;
    int neles=0;
    
    RefinedGeometry *refined_geo;
    Array<int>       coloring;
    
    //
    // find the # of output points and cells at the selected level of 
    // refinement (we may want to cache this...)
    //
    for (int i = 0; i < mesh->GetNE(); i++)
    {
        int geom = mesh->GetElementBaseGeometry(i);
        int ele_nverts = Geometries.GetVertices(geom)->GetNPoints();
        refined_geo    = GlobGeometryRefiner.Refine(geom, lod, 1);
        npts  += refined_geo->RefPts.GetNPoints();
        neles += refined_geo->RefGeoms.Size() / ele_nverts;
    }

    vtkFloatArray *rv = vtkFloatArray::New();
    rv->SetNumberOfComponents(1);
    rv->SetNumberOfTuples(neles);

    // set output array value from the mfem mesh's "Attribue" field
    int ref_idx=0;
    for (int i = 0; i < mesh->GetNE(); i++)
    {
        int geom = mesh->GetElementBaseGeometry(i);
        int nv = Geometries.GetVertices(geom)->GetNPoints();
        refined_geo= GlobGeometryRefiner.Refine(geom, lod, 1);
        int attr = mesh->GetAttribute(i);
        for (int j = 0; j < refined_geo->RefGeoms.Size(); j += nv)
        {
             rv->SetTuple1(ref_idx,attr);
             ref_idx++;
        }
   }
   
   delete mesh;
   
   return rv;
}


// ****************************************************************************
//  Method: avtMFEMFileFormat::RegisterDataSelections
//
//  Purpose: 
//   Used to suport avtResolutionSelection & capture the selected lod.
//
//  Arguments:
//     sels:    data selection list from the pipeline
//     applied: pipeline handshaking for handling data selections
//
//
// ****************************************************************************
void
avtMFEMFileFormat::RegisterDataSelections(const std::vector<avtDataSelection_p>& sels, 
                                          std::vector<bool>* applied)
{
    for(size_t i=0; i < sels.size(); ++i)
    {
        if(strcmp(sels[i]->GetType(), "avtResolutionSelection") == 0)
        {
            const avtResolutionSelection* sel =
                static_cast<const avtResolutionSelection*>(*sels[i]);
            this->selectedLOD = sel->resolution();
            (*applied)[i] = true;
        }
    }
}



