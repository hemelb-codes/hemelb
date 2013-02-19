//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELBSETUPTOOL_POLYDATAGENERATOR_H
#define HEMELBSETUPTOOL_POLYDATAGENERATOR_H

#include "GeometryGenerator.h"

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

class PolyDataGenerator: public GeometryGenerator {
public:
	PolyDataGenerator();
	virtual ~PolyDataGenerator();

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
	virtual void ComputeBounds(double[]) const;
	virtual void PreExecute(void);
	void ClassifySite(Site& site);

	/**
	 * This method implements the algorithm used to approximate the wall normal at a given
	 * fluid site. This is done based on the normal of the triangles intersected by
	 * each lattice link and the distance to those intersections. The method will return
	 * whether it was possible to compute such normal (i.e. any intersecting link?).
	 *
	 * Current implementation does a weighted sum of the wall normals. The weights are the
	 * reciprocal of cut distances along each link.
	 *
	 * @param normal Computed normal.
	 * @param intersectedCellsNormals Vector with the normals of all the triangles intersected by a link.
	 * @param intersectedCellsDistance Vector with the distances to the intersections along a link.
	 * @return Whether it was possible to compute a wall normal.
	 */
	bool ComputeAveragedNormal(Vector& normal,
			std::vector<double*>& intersectedCellsNormals,
			std::vector<float>& intersectedCellsDistance) const;

	int ComputeIntersections(Site& from, Site& to);
	// represents whether the block is inside (-1) outside (+1) or undetermined (0)
	virtual int BlockInsideOrOutsideSurface(const Block &block);
	// Members set from outside to initialise
	double SeedPointWorking[3];
	vtkPolyData* ClippedSurface;
	vtkOBBTree* Locator;

	// Members used internally
	vtkPoints* hitPoints;
	vtkIdList* hitCellIds;
	vtkIntArray* IoletIdArray;
};

#endif // HEMELBSETUPTOOL_POLYDATAGENERATOR_H
