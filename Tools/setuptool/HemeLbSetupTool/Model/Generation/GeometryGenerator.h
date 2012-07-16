// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELBSETUPTOOL_GEOMETRYGENERATOR_H
#define HEMELBSETUPTOOL_GEOMETRYGENERATOR_H

#include <string>
#include <vector>

// VTK bits we need
class vtkPolyData;
class vtkOBBTree;
class vtkPoints;
class vtkIdList;
class vtkIntArray;

#include "GetSet.h"
#include "Iolet.h"
#include "GenerationError.h"

class GeometryWriter;
class Site;
class BlockWriter;

struct InconsistentFluidnessError: public GenerationError {
	InconsistentFluidnessError(const Site& s1, const Site& s2, const int nHits);
	~InconsistentFluidnessError() throw () {
	}
	virtual const char* what() const throw ();
	const Site& s1;
	const Site& s2;
	const int nHits;
private:
	std::string msg;
};

class GeometryGenerator {
public:
	GeometryGenerator();
	~GeometryGenerator();
	void Execute() throw (GenerationError);

	inline double GetVoxelSizeMetres(void) {
		return this->VoxelSizeMetres;
	}
	inline void SetVoxelSizeMetres(double val) {
		this->VoxelSizeMetres = val;
	}

	inline double GetVoxelSizeWorking(void) {
		return 1;
	}

	inline std::string GetOutputGeometryFile(void) {
		return this->OutputGeometryFile;
	}
	inline void SetOutputGeometryFile(std::string val) {
		this->OutputGeometryFile = val;
	}

	inline std::vector<Iolet*>& GetIolets() {
		return this->Iolets;
	}
	inline void SetIolets(std::vector<Iolet*> iv) {
		this->Iolets = std::vector<Iolet*>(iv);
	}

	inline void GetSeedPointWorking(double out[3]) {
		for (unsigned int i = 0; i < 3; ++i)
			out[i] = this->SeedPointWorking[i];
		return;
	}
	inline void SetSeedPointWorking(double out[3]) {
		for (unsigned int i = 0; i < 3; ++i)
			this->SeedPointWorking[i] = out[i];
	}
	inline void SetSeedPointWorking(double x, double y, double z) {
		this->SeedPointWorking[0] = x;
		this->SeedPointWorking[1] = y;
		this->SeedPointWorking[2] = z;
	}

	inline vtkPolyData* GetClippedSurface(void) {
		return this->ClippedSurface;
	}
	inline void SetClippedSurface(vtkPolyData* val) {
		this->ClippedSurface = val;
	}

private:
	void ClassifySite(Site& site);
	void WriteSolidSite(BlockWriter& blockWriter, Site& site);
	void WriteFluidSite(BlockWriter& blockWriter, Site& site);
	int ComputeIntersections(Site& from, Site& to);
	// Members set from outside to initialise
	double VoxelSizeMetres;
	std::string OutputGeometryFile;
	std::vector<Iolet*> Iolets;
	double SeedPointWorking[3];
	vtkPolyData* ClippedSurface;
	vtkOBBTree* Locator;

	// Members used internally
	vtkPoints* hitPoints;
	vtkIdList* hitCellIds;
	vtkIntArray* IoletIdArray;
};

#endif // HEMELBSETUPTOOL_GEOMETRYGENERATOR_H
