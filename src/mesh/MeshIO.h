#pragma once

#include <MaterialTypes.h>

namespace g3 {


enum class IOCode
{
	Ok = 0,

	FileAccessError = 1,
	UnknownFormatError = 2,
	FormatNotSupportedError = 3,
	InvalidFilenameError = 4,

	// read errors
	FileParsingError = 100,
	GarbageDataError = 101,
	GenericReaderError = 102,
	GenericReaderWarning = 103,

	// write errors
	WriterError = 200,


	// other status
	ComputingInWorkerThread = 1000
};



class ReadOptions
{
public:
	bool ReadMaterials;

	// [RMSC#] disabled
	// format readers will inevitably have their own settings, we
	// can use this to pass arguments to them
	//CommandArgumentSet CustomFlags = new CommandArgumentSet();

	ReadOptions()
	{
		ReadMaterials = false;
	}

	static ReadOptions Defaults() {
		ReadOptions r;
		r.ReadMaterials = false;
		return r;
	};
};




struct IOReadResult
{
	IOCode code;
	std::string message;

	IOReadResult(IOCode r, const std::string & s)
	{
		code = r;
		message = s;
		if (message == "")
			message = "(no message)";
	}

	static IOReadResult Ok() {
		return IOReadResult(IOCode::Ok, "");
	}
};







struct IOWriteResult
{
	IOCode code;
	std::string message;

	IOWriteResult(IOCode r, const std::string & s)
	{
		code = r;
		message = s;
		if (message == "")
			message = "(no message)";
	}

	static IOWriteResult Ok() {
		return IOWriteResult(IOCode::Ok, "");
	}
};




struct WriteOptions
{
	bool bWriteBinary;        	// write binary format if supported (STL)

	bool bPerVertexNormals;		// write per-vertex normals (OBJ)
	bool bPerVertexColors;		// write per-vertex colors (OBJ)
	bool bPerVertexUVs;			// write per-vertex UVs
										// can be overridden by per-mesh UVs in WriteMesh
	bool bWriteGroups;			// write face groups (OBJ)

	bool bCombineMeshes;     	// combine all input meshes into a single output mesh
										// some STL readers do not handle multiple solids...

	int RealPrecisionDigits;		// number of digits of float precision (after decimal)

	bool bWriteMaterials;		// for OBJ, indicates that .mtl file should be written
	std::string MaterialFilePath;		// only used if bWriteMaterialFile = true

	std::string groupNamePrefix;        // prefix for group names in OBJ files (default is "mmGroup")
	std::function<std::string(int)> GroupNameF;  // if non-null, you can use this to generate your own group names


	std::function<void(int,int)> ProgressFunc;	// progress monitoring callback

	std::function<std::string()> AsciiHeaderFunc;    // if you define this, returned string will be written as header start of ascii formats

	static WriteOptions Defaults() {
		WriteOptions opt;
		opt.bWriteBinary = false;
		opt.bPerVertexNormals = false;
		opt.bPerVertexColors = false;
		opt.bWriteGroups = false;
		opt.bPerVertexUVs = false;
		opt.bCombineMeshes = false;
		opt.bWriteMaterials = false;
		opt.ProgressFunc = nullptr;
		opt.RealPrecisionDigits = 15;
		//opt.RealPrecisionDigits = 7;
		return opt;
	};
};




class DenseUVMesh
{
public:
	dvector<Vector2f> UVs;
	dvector<Index3i> TriangleUVs;
};



struct WriteMesh
{
	DMesh3Ptr Mesh;
	std::string Name;         // supported by some formats

	// [RMS] disable this for now

	std::vector<GenericMaterialPtr> Materials;		// set of materials (possibly) used in this mesh
	std::map<int, int> TriToMaterialMap;		// triangle index -> Materials list index
	//IIndexMap TriToMaterialMap;			

	DenseUVMesh * UVs;  // separate UV layer (just one for now)
							 // assumption is that # of triangles in this UV mesh is same as in Mesh

	WriteMesh(DMesh3Ptr mesh, std::string name = "") : Mesh(mesh) {
		Name = name;
		//UVs = nullptr;
		//Materials = nullptr;
		//TriToMaterialMap = nullptr;
	}
};





}