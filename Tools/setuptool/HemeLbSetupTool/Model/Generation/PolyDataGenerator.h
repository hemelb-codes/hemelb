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

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Verbose_ostream.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 PointCGAL;
typedef Kernel::Plane_3 PlaneCGAL;
typedef Kernel::Vector_3 VectorCGAL;
typedef Kernel::Segment_3 SegmentCGAL;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Polyhedron::Vertex_iterator     Vertex_iteratorCGAL;
typedef Tree::Object_and_primitive_id Object_and_primitive_id;
typedef Tree::Primitive_id Primitive_id;


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
	inline void SetClippedCGALSurface(Polyhedron* val) {
		this->ClippedCGALSurface = val;
	}

private:
	virtual void ComputeBounds(double[]) const;
	virtual void PreExecute(void);
	void ClassifySite(Site& site);
	int ComputeIntersections(Site& from, Site& to);
	// represents whether the block is inside (-1) outside (+1) or undetermined (0)
	virtual int BlockInsideOrOutsideSurface(const Block &block);
	// Members set from outside to initialise
	double SeedPointWorking[3];
	vtkPolyData* ClippedSurface;
	Polyhedron* ClippedCGALSurface;
	vtkOBBTree* Locator;
	Tree* AABBtree;

	// Members used internally
	vtkPoints* hitPoints;
	vtkIdList* hitCellIds;
	vtkIntArray* IoletIdArray;
};

#endif // HEMELBSETUPTOOL_POLYDATAGENERATOR_H
