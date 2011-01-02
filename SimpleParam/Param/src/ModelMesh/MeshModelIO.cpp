/* ================== Library Information ================== */
// [Name] 
// MeshLib Library
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// A general, flexible, versatile and easy-to-use mesh library for research purpose.
// Supporting arbitrary polygonal meshes as input, but with a 
// primary focus on triangle meshes.

/* ================== File Information ================== */
// [Name]
// MeshModelIO.cpp
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// Defining the input/output (I/O) functions of the kernel components of the mesh model
// Supporting various 3D model files, including .tm, .obj, .off, etc

#include "MeshModelIO.h"
#include <iostream>
#include <fstream>

/* ================== Mesh Model I/O Functions ================== */

// Constructor
MeshModelIO::MeshModelIO()
{
    kernel = NULL;
}

// Destructor
MeshModelIO::~MeshModelIO()
{
}

// Initializer
void MeshModelIO::ClearData()
{

}

void MeshModelIO::AttachKernel(MeshModelKernel* pKernel /* = NULL */)
{
    kernel = pKernel;
}



// General I/O functions, platform-dependent functions 
bool MeshModelIO::LoadModel(const std::string& filename)
{
    // Resolve file name
    std::string file_path, file_title, file_ext;
    util.ResolveFileName(filename, file_path, file_title, file_ext);

    printf("Load Model %s... ", (file_title+file_ext).c_str());

    // Set model informaion
    ModelInfo& mInfo = kernel->GetModelInfo();
    mInfo.SetFileName(filename);

    // Load model
    util.MakeLower(file_ext);

    bool bOpenFlag;
    if(file_ext == ".tm")
    {
        bOpenFlag = OpenTmFile(filename);
    }
    else if(file_ext == ".ply2")
    {
        bOpenFlag = OpenPly2File(filename);
    }
    else if(file_ext == ".off")
    {
        bOpenFlag = OpenOffFile(filename);
    }
	else if(file_ext == ".obj")
    {
        bOpenFlag = OpenObjFile(filename);
    }
    else
    {
        bOpenFlag = false;
    }

    printf("%s\n", (bOpenFlag == true) ? "Success" : "Fail");
    if(bOpenFlag)
    {
        size_t nVertex = kernel->GetVertexInfo().GetCoord().size();
        size_t nFace   = kernel->GetFaceInfo().GetIndex().size();
        printf("#Vertex = %d, #Face = %d\n\n", nVertex, nFace);
    }

    return bOpenFlag;
}

bool MeshModelIO::StoreModel(const std::string& filename)
{
    // Resolve file name
    std::string file_path, file_title, file_ext;
    util.ResolveFileName(filename, file_path, file_title, file_ext);

    printf("Store Model %s... ", (file_title+file_ext).c_str());

    // Store model
    util.MakeLower(file_ext);

    bool bSaveFlag;
    if(file_ext == ".tm")
    {
        bSaveFlag = SaveTmFile(filename);
    }
    else if(file_ext == ".ply2")
    {
        bSaveFlag = SavePly2File(filename);
    }
	else if(file_ext == ".off")
    {
        bSaveFlag = SaveOffFile(filename);
    }
	else if(file_ext == ".obj")
    {
        bSaveFlag = SaveObjFile(filename);
    }
    else
    {
        bSaveFlag = false;
    }

    printf("%s\n", (bSaveFlag == true) ? "Success" : "Fail");

    return bSaveFlag;
}

// .tm file I/O functions
bool MeshModelIO::OpenTmFile(const std::string& filename)
{
    // Read data from file
    FILE* fp = fopen(filename.c_str(), "r");    

    if(fp == NULL)
        return false;

    int nVertex, nFace, nReadFace;
    fscanf(fp, "%d", &nVertex);
    fscanf(fp, "%d", &nFace);
    fscanf(fp, "%d", &nReadFace);

    // Prepare for a new mesh model
    kernel->ClearData();

    // Load vertex information
    VertexInfo& vInfo = kernel->GetVertexInfo();
    CoordArray& vCoord = vInfo.GetCoord();

    vCoord.resize(nVertex);

    int i;
    float coordPos0;
    float coordPos1;
    float coordPos2;
    for(i = 0; i < nVertex; ++ i)
    {
        Coord& v = vCoord[i];
        fscanf(fp, "%f %f %f", &coordPos0, &coordPos1, &coordPos2);
        v[0] = coordPos0;
        v[1] = coordPos1;
        v[2] = coordPos2;
    }

    // Load face information
    FaceInfo& fInfo = kernel->GetFaceInfo();
    PolyIndexArray& fIndex = fInfo.GetIndex();

    fIndex.resize(nFace);
    int fv0;
    int fv1;
    int fv2;
    for(i = 0; i < nFace; ++ i)
    {
        IntArray& f = fIndex[i];
        f.resize(3);

        fscanf(fp, "%d %d %d", &fv0, &fv1, &fv2);
        f[0] = fv0;
        f[1] = fv1;
        f[2] = fv2;
    }

    fclose(fp);

    return true;
}

bool MeshModelIO::SaveTmFile(const std::string& filename)
{
    // Read data from file
    std::ofstream file(filename.c_str());
    if(!file)
        return false;

    size_t nVertex, nFace;
    VertexInfo& vInfo = kernel->GetVertexInfo();
    CoordArray& arrCoord = vInfo.GetCoord();
    nVertex = arrCoord.size();

    FaceInfo& fInfo = kernel->GetFaceInfo();
    PolyIndexArray& arrIndex = fInfo.GetIndex();
    nFace = arrIndex.size();

    file << nVertex << ' ' << nFace << ' ' << nFace << '\n';

    // Store vertex information
    size_t i, j;
    for(i = 0; i < nVertex; ++ i)
    {
        Coord& v = arrCoord[i];
        for(j = 0; j < 3; ++ j)
            file << v[j] << ((j<2) ? ' ' : '\n');
    }

    // Store face information
    for(i = 0; i < nFace; ++ i)
    {
        IntArray& f = arrIndex[i];
        for(j = 0; j < 3; ++ j)
            file << f[j] << ((j<2) ? ' ' : '\n');
    }
    file.close();

    return true;
}

// .ply2 file I/O functions
bool MeshModelIO::OpenPly2File(const std::string& filename)
{
    // Read data from file
    std::ifstream file(filename.c_str());
    if(! file)
        return false;

    int nVertex, nFace;
    file >> nVertex >> nFace;

    // Prepare for a new mesh model
    kernel->ClearData();

    // Load vertex information
    VertexInfo& vInfo = kernel->GetVertexInfo();
    CoordArray& vCoord = vInfo.GetCoord();

    vCoord.resize(nVertex);
    int i, j, n;
    for(i = 0; i < nVertex; ++ i)
    {
        Coord& v = vCoord[i];
        for(j = 0; j < 3; ++ j)
            file >> v[j];
    }

    // Load face information
    FaceInfo& fInfo = kernel->GetFaceInfo();
    PolyIndexArray& fIndex = fInfo.GetIndex();

    fIndex.resize(nFace);
    for(i = 0; i < nFace; ++ i)
    {
        IntArray& f = fIndex[i];
        file >> n;
        f.resize(n);
        for(j = 0; j < n; ++ j)
            file >> f[j];
    }
    file.close();

    return true;
}

bool MeshModelIO::SavePly2File(const std::string& filename)
{
    // Read data from file
    std::ofstream file(filename.c_str());
    if(!file)
        return false;

    size_t nVertex, nFace;
    VertexInfo& vInfo = kernel->GetVertexInfo();
    CoordArray& arrCoord = vInfo.GetCoord();
    nVertex = arrCoord.size();

    FaceInfo& fInfo = kernel->GetFaceInfo();
    PolyIndexArray& arrIndex = fInfo.GetIndex();
    nFace = arrIndex.size();

    file << nVertex << ' ' << nFace << '\n';

    // Store vertex information
    size_t i, j, n;
    for(i = 0; i < nVertex; ++ i)
    {
        Coord& v = arrCoord[i];
        for(j = 0; j < 3; ++ j)
            file << v[j] << ((j<2) ? ' ' : '\n');
    }

    // Store face information
    for(i = 0; i < nFace; ++ i)
    {
        IntArray& f = arrIndex[i];
        n = f.size();
        file << n << ' ';
        for(j = 0; j < 3; ++ j)
            file << f[j] << ((j<2) ? ' ' : '\n');
    }
    file.close();

    return true;
}

// .off file I/O functions
bool MeshModelIO::OpenOffFile(const std::string& filename)
{
    // Read data from file
    FILE* fp = fopen(filename.c_str(), "r");

    if(fp == NULL)
        return false;
    char format[10];
    fscanf(fp, "%s", format);

    int nVertex, nFace, zero;
    fscanf(fp, "%d", &nVertex);
    fscanf(fp, "%d", &nFace);
    fscanf(fp, "%d", &zero);

    // Prepare for a new mesh model
    kernel->ClearData();

    // Load vertex information
    VertexInfo& vInfo = kernel->GetVertexInfo();
    CoordArray& vCoord = vInfo.GetCoord();

    vCoord.resize(nVertex);
    int i, j, n;
    for(i = 0; i < nVertex; ++ i)
    {
        float tempCoord;
        Coord& v = vCoord[i];
        for(j = 0; j < 3; ++ j)
        {
            fscanf(fp, "%f", &tempCoord);
            v[j] = tempCoord;
        }
    }

    // Load face information
    FaceInfo& fInfo = kernel->GetFaceInfo();
    PolyIndexArray& fIndex = fInfo.GetIndex();

    fIndex.resize(nFace);
    for(i = 0; i < nFace; ++ i)
    {
        IntArray& f = fIndex[i];
        fscanf(fp, "%d", &n);
        f.resize(n);

        int vertexID;
        for(j = 0; j < n; ++ j)
        {
            fscanf(fp, "%d", &vertexID);
            f[j] = vertexID;
        }
    }
    fclose(fp);

    return true;
}

// .off file I/O functions
bool MeshModelIO::SaveOffFile(const std::string& filename)
{
    // Read data from file
    std::ofstream file(filename.c_str());
    if(!file)
        return false;

	file << "OFF" << std::endl;
	
    size_t nVertex, nFace;
    VertexInfo& vInfo = kernel->GetVertexInfo();
    CoordArray& arrCoord = vInfo.GetCoord();
    nVertex = arrCoord.size();
    FaceInfo& fInfo = kernel->GetFaceInfo();
    PolyIndexArray& arrIndex = fInfo.GetIndex();
    nFace = arrIndex.size();
	
    file << nVertex << ' ' << nFace << ' ' << 0 <<'\n';
	
	// Store vertex information
    size_t i, j, n;
    for(i = 0; i < nVertex; ++ i)
    {
        Coord& v = arrCoord[i];
        for(j = 0; j < 3; ++ j)
            file << v[j] << ((j<2) ? ' ' : '\n');
    }
	
    // Store face information
    for(i = 0; i < nFace; ++ i)
    {
        IntArray& f = arrIndex[i];
        n = f.size();
        file << n << ' ';
        for(j = 0; j < 3; ++ j)
            file << f[j] << ((j<2) ? ' ' : '\n');
    }
    file.close();
	
    return true;
}
bool MeshModelIO::OpenObjFile(const std::string& filename)
{
	int nVertex = 0, nFace = 0;
	float coordPos0;
    float coordPos1;
    float coordPos2;
	int fv0;
    int fv1;
    int fv2;
	char format[10];


    std::string str;
    std::ifstream ifs(filename.c_str());

	if(ifs.fail()) return false;

	// find the vertex and face number here.
	while(!ifs.eof()) {
		getline(ifs, str, '\n');
		if(ifs.fail()) break;
		if(str.empty()) continue;
		

		if(str[0] == '#')	continue;
		if(str[0] == 'g')	continue;
		if(str[0] == 'v' && str[1] == ' ') 
		{
			++nVertex;
			continue;
		}
		if(str[0] == 'f' && str[1] == ' ') 
		{
			++nFace;
			continue;
		}
	}

	ifs.clear();
	ifs.seekg(0, std::ifstream::beg);

	// Prepare for a new mesh model
    kernel->ClearData();
	
    // Load vertex information
    VertexInfo& vInfo = kernel->GetVertexInfo();
    CoordArray& vCoord = vInfo.GetCoord();
	
    vCoord.resize(nVertex);

	// Load Texture information
	TexCoordArray&  tCoord = vInfo.GetTexCoord();
	tCoord.resize(nVertex);


	// Load face information
    FaceInfo& fInfo = kernel->GetFaceInfo();
    PolyIndexArray& fIndex = fInfo.GetIndex();

    fIndex.resize(nFace);


	size_t vn = 0, fn = 0, vt=0;

	while(!ifs.eof()) {
		getline(ifs, str, '\n');
		if(ifs.fail()) break;
		if(str.empty()) continue;
		
		
		if(str[0] == '#')	continue;
		if(str[0] == 'g')	continue;
		if(str[0] == 'o')	continue;

		// read the vertex information here.
		if(str[0] == 'v' && str[1] == ' ') 
		{
			sscanf(str.c_str(), "%s %f %f %f", format, &coordPos0, &coordPos1, &coordPos2);

			Coord& v = vCoord[vn];
			v[0] = coordPos0;
			v[1] = coordPos1;
            v[2] = coordPos2;

			++vn;
			continue;
		}

		if(str[0] == 'v' && str[1] == 't')
		{
			sscanf(str.c_str(), "%s %f %f", format, &coordPos0, &coordPos1);

			TexCoord& v = tCoord[vt];
			v[0] = coordPos0;
			v[1] = coordPos1;
			
			++vt;
			continue;
		}

		// read the face information here.
		if(str[0] == 'f' && str[1] == ' ') 
		{
            std::vector<int> word_pos(1, 0);
            std::string delimit = " \t";
			for(size_t i = 0; i < str.length(); ++i) {
				if(delimit.find(str[i]) == std::string::npos)
					continue;
				if(i+1 < str.length() && delimit.find(str[i+1]) == std::string::npos)
					word_pos.push_back((int) i+1);
			}
			if(word_pos.size() < 4) {
                std::cout << "bad polygon";
				return false;
			}
			sscanf(str.c_str()+word_pos[1], "%d", &fv0);
			sscanf(str.c_str()+word_pos[2], "%d", &fv1);
			sscanf(str.c_str()+word_pos[3], "%d", &fv2);
						
			IntArray& f = fIndex[fn];
			f.resize(3);
			
			f[0] = fv0 - 1;
			f[1] = fv1 - 1;
            f[2] = fv2 - 1;

			++fn;

			continue;
		}

	}

	ifs.close();

	return true;
}
bool MeshModelIO::SaveObjFile(const std::string& filename)
{
    std::ofstream file(filename.c_str());
    if(!file)
        return false;

	file << "# " << std::endl;
	file << "# Wavefront OBJ file" << std::endl;
	file << "# object ..." + filename << std::endl;
	
    size_t nVertex, nFace;
    VertexInfo& vInfo = kernel->GetVertexInfo();
    CoordArray& arrCoord = vInfo.GetCoord();
    nVertex = arrCoord.size();
	
    FaceInfo& fInfo = kernel->GetFaceInfo();
    PolyIndexArray& arrIndex = fInfo.GetIndex();
    nFace = arrIndex.size();
	
	// Store vertex information
    size_t i, j;
    for(i = 0; i < nVertex; ++ i)
    {
        Coord& v = arrCoord[i];
		file << "v ";
        for(j = 0; j < 3; ++ j)
            file << v[j] << ((j<2) ? ' ' : '\n');
    }
	
    // Store face information
    for(i = 0; i < nFace; ++ i)
    {
        IntArray& f = arrIndex[i];
        file << "f ";
        for(j = 0; j < 3; ++ j)
            file << f[j] + 1 << ((j<2) ? ' ' : '\n');
    }

    file.close();
	
	return true;
}


/* \function LoadTriangularMesh
 * load the triangular mesh from a coord-array and face-list
 */
bool MeshModelIO::LoadTriangularMesh(boost::shared_ptr<const tri_mesh_3d> p_mesh)
{
	if(p_mesh == NULL) 
		return false;

	size_t nVertex = (size_t) (p_mesh->vert_num);
	size_t nFace = (size_t) (p_mesh->face_num);
	
	double* vertex_coord_array = p_mesh->vertex;
	int* face_list_array = p_mesh->face;

	double x_coord, y_coord, z_coord;

	// Prepare for a new mesh model
	kernel->ClearData();

	// Load vertex information
	VertexInfo& vInfo = kernel->GetVertexInfo();
	CoordArray& vCoord = vInfo.GetCoord();

	vCoord.resize(nVertex);
	for(size_t vn = 0; vn < nVertex; ++vn)
	{
		x_coord = vertex_coord_array[vn*3 + 0];
		y_coord = vertex_coord_array[vn*3 + 1];
		z_coord = vertex_coord_array[vn*3 + 2];
		vCoord[vn] = Coord(x_coord, y_coord, z_coord);
	}

	// Load face information
	FaceInfo& fInfo = kernel->GetFaceInfo();
	PolyIndexArray& fIndex = fInfo.GetIndex();

	fIndex.resize(nFace);
	for(size_t fn = 0; fn < nFace; ++fn)
	{
		IntArray& f = fIndex[fn];
		for(size_t k=0; k<3; ++k)
		{
			f[k] = face_list_array[fn*3+k];
		}
	}

	return true;
}
