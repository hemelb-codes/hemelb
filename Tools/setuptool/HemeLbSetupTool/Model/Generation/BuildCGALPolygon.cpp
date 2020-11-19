// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "BuildCGALPolygon.h"
#include "io/formats/geometry.h"


#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

using namespace hemelb::io::formats;

template<class HDS> 
void BuildCGALPolygon<HDS>::operator()( HDS& hds){
	CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
	B.begin_surface( this->pts->GetNumberOfPoints(), this->polys->GetNumberOfCells(), 0);
	typedef typename HDS::Vertex   Vertex;
	typedef typename Vertex::Point Point;
	vtkIdType npts = 0;
	vtkIdType *indx;
	double vertex[3];
	for (int i = 0; i != this->pts->GetNumberOfPoints(); i++ ){
		this->pts->GetPoint(i, vertex);
		B.add_vertex( Point( vertex[0], vertex[1], vertex[2]));

	}
	int j = 0;
	Face_handle face; 
	for (this->polys->InitTraversal(); this->polys->GetNextCell(npts,indx); ){
		if(indx[0] != indx[1] & indx[0] != indx[2] & indx[1] != indx[2]){ 
			//VTK polygons can contain lines where two vertexes are identical. Forget these
			if (B.test_facet(indx, indx+3)){ 
				face = B.begin_facet();
				B.add_vertex_to_facet( indx[0]);
				B.add_vertex_to_facet( indx[1]);
				B.add_vertex_to_facet( indx[2]);
				B.end_facet();
				//cout << this->IoletIdArray->GetValue(j) << endl;
				face->id() = IoletIdArray->GetValue(j) + 2;
				//the face id is size_t i.e. unsigned so we shift this to positive. 1 is wall. 2,3 ... are 
				//the inlets and outlets.  
			}
			else
				cout << "Ignoring Non manifold facet between: " << indx[0] << " " <<  indx[1] << " " << indx[2] << endl;
		}
		else{
			cout << "Eleminated degenerate vertex: " << indx[0] << " " <<  indx[1] << " " << indx[2] << endl;
		}
		++j;
	}
	
	B.end_surface();
	//cout << B.check_unconnected_vertices () << endl;
	B.remove_unconnected_vertices();
	//cout << B.check_unconnected_vertices () << endl;
}

