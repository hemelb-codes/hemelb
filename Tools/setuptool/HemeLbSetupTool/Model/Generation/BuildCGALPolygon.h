//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELBSETUPTOOL_BUILDCGALPOLYGON_H
#define HEMELBSETUPTOOL_BUILDCGALPOLYGON_H


#include "CGALtypedef.h"

#include "vtkIdList.h"
#include "vtkCellData.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkOBBTree.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
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
	BuildCGALPolygon(vtkPoints* ptsin, vtkCellArray *polysin) {
		this->pts = ptsin;
		this->polys = polysin;
	}
	void operator()( HDS& hds);

private:
	vtkPoints* pts;
	vtkCellArray* polys;
	
};


template<class HDS> void BuildCGALPolygon<HDS>::operator()( HDS& hds){
	CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
	B.begin_surface( this->pts->GetNumberOfPoints(), this->polys->GetNumberOfCells(), 0);
	typedef typename HDS::Vertex   Vertex;
	typedef typename Vertex::Point Point;
	vtkIdType npts = 0;
	vtkIdType *indx = 0;
	double vertex[3];
	for (int i = 0; i != this->pts->GetNumberOfPoints(); i++ ){
		this->pts->GetPoint(i, vertex);
		//cout << i << " " << vertex[0] << " " <<vertex[1] << " " << vertex[2] << endl;
		B.add_vertex( Point( vertex[0], vertex[1], vertex[2]));

	}
	int k = 0;
	for (this->polys->InitTraversal(); this->polys->GetNextCell(npts,indx); ){
		if(indx[0] != indx[1] & indx[0] != indx[2] & indx[1] != indx[2]){ 
			//VTK polygons can contain lines where two vertexes are identical. Forget these
			B.begin_facet();
			B.add_vertex_to_facet( indx[0]);
			B.add_vertex_to_facet( indx[1]);
			B.add_vertex_to_facet( indx[2]);
			//cout << k << " " << indx[0] << " "<< indx[1] << " " << indx[2] << " " << endl;
			B.end_facet();
			++k;
		}
		else{
			//We need to acout for this in the iolet map.
			cout << "Eleminated degenrate vertex" << endl;
		}
	}	
	B.end_surface();
}

#endif //HEMELBSETUPTOOL_BUILDCGALPOLYGON_H
