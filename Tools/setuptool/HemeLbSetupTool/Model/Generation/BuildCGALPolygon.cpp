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


template<class HDS> void BuildCGALPolygon<HDS>::operator()( HDS& hds){
	CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
	B.begin_surface( this->pts->GetNumberOfPoints(), this->polys->GetNumberOfCells(), 0);
	typedef typename HDS::Vertex   Vertex;
	typedef typename Vertex::Point Point;
	vtkIdType npts = 0;
	vtkIdType *indx;
	double vertex[3];
	for (int i = 0; i != this->pts->GetNumberOfPoints(); i++ ){
		this->pts->GetPoint(i, vertex);
		//cout << i << " " << vertex[0] << " " <<vertex[1] << " " << vertex[2] << endl;
		B.add_vertex( Point( vertex[0], vertex[1], vertex[2]));

	}
	//int k = 0;
	int j = 0;
	for (this->polys->InitTraversal(); this->polys->GetNextCell(npts,indx); ){
		if(indx[0] != indx[1] & indx[0] != indx[2] & indx[1] != indx[2]){ 
			//VTK polygons can contain lines where two vertexes are identical. Forget these
			std::vector<vtkIdType>  testvec;
			testvec.push_back(indx[0]);
			testvec.push_back(indx[1]);
			testvec.push_back(indx[2]);
			if (B.test_facet(testvec.begin(), testvec.end())){
				B.begin_facet();
				B.add_vertex_to_facet( indx[0]);
				B.add_vertex_to_facet( indx[1]);
				B.add_vertex_to_facet( indx[2]);
				//cout << k << " " << indx[0] << " "<< indx[1] << " " << indx[2] << " " << endl;
				B.end_facet();
				//cout << this->IoletIdArray->GetValue(j) << endl;
				ID.push_back(this->IoletIdArray->GetValue(j));
			}
			else
				cout << "ignoring " << testvec[0] << " " <<  testvec[1] << " " << testvec[2] << endl;
		}
		else{
			//We need to acout for this in the iolet map.
			cout << "Eleminated degenerate vertex" << endl;
		}
		++j;
	}	
	B.end_surface();
}
