// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELBSETUPTOOL_BUILDCGALPOLYGON_H
#define HEMELBSETUPTOOL_BUILDCGALPOLYGON_H


#include "CGALtypedef.h"

#include "vtkIdList.h"
#include "vtkCellData.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkOBBTree.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"
#include "vtkMatrix4x4.h"
#include "Block.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkCellArray.h"


#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>



template <class HDS>
class BuildCGALPolygon: public CGAL::Modifier_base<HDS> {
public:
	BuildCGALPolygon(vtkPoints* ptsin, vtkCellArray *polysin,vtkIntArray* IoletIdArrayIn) {
		this->pts = ptsin;
		this->polys = polysin;
		this->IoletIdArray = IoletIdArrayIn;
	}
	void operator()( HDS& hds);
	
private:
	vtkPoints* pts;
	vtkCellArray* polys;
	vtkIntArray* IoletIdArray;
	
};

#endif //HEMELBSETUPTOOL_BUILDCGALPOLYGON_H
