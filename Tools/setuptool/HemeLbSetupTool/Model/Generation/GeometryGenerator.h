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

#include "Iolet.h"
#include "GenerationError.h"

class GeometryWriter;
class Site;
class BlockWriter;
class Block;

class GeometryGenerator {
public:
	GeometryGenerator();
	virtual ~GeometryGenerator();
	void Execute(bool skipNonIntersectingBlocks) throw (GenerationError);

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

protected:
	virtual void ComputeBounds(double []) const = 0;
	virtual void PreExecute(void);
	virtual void ClassifySite(Site& site) = 0;
	void WriteSolidSite(BlockWriter& blockWriter, Site& site);
	void WriteFluidSite(BlockWriter& blockWriter, Site& site);
	// Members set from outside to initialise
	double VoxelSizeMetres;
	std::string OutputGeometryFile;
	std::vector<Iolet*> Iolets;
    virtual bool BlockIntersectsSurface(const Block &block,int & side) = 0;
};

#endif // HEMELBSETUPTOOL_GEOMETRYGENERATOR_H
