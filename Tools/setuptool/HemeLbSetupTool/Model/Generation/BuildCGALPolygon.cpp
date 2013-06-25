//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "BuildCGALPolygon.h"
#include "io/formats/geometry.h"


#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

using namespace hemelb::io::formats;


BuildCGALPolygon::BuildCGALPolygon(vtkPoints* ptsin, vtkCellArray *polysin) {
		this->pts = ptsin;
		this->polys = polysin;
	}

void BuildCGALPolygon::testing( HDS& hds) {
// Postcondition: hds is a valid polyhedral surface.
CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
	B.begin_surface( 86, 168, 0);
	typedef typename HDS::Vertex   Vertex;
	typedef typename Vertex::Point Point;
	vtkIdType npts = 0;
	vtkIdType *indx = 0;
	double n[3], v1[3], v2[3], v3[3];
	for (int i = 0; i != 86; i++ ){
		this->pts->GetPoint(i,v1);
		this->pts->GetPoint(i,v2);
		this->pts->GetPoint(i,v3);
		cout << i << " " << v1[0] << " " <<v1[1] << " " << v1[2] << endl;
B.add_vertex( Point( v1[0], v1[1], v1[2]));
B.add_vertex( Point( v2[0], v2[1], v2[2]));
		B.add_vertex( Point( v3[0], v3[1], v3[2]));
	}
	int k = 0;
		for (this->polys->InitTraversal(); this->polys->GetNextCell(npts,indx); ){

			B.begin_facet();
			B.add_vertex_to_facet( indx[0]);
			B.add_vertex_to_facet( indx[1]);
			B.add_vertex_to_facet( indx[2]);
			cout << k << " " << indx[0] << " "<< indx[1] << " " << indx[2] << " " << endl;
			B.end_facet();
			++k;
		}
		B.end_surface();
}




