//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#include "MkCGALMesh.h"
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

class SurfaceCreator: public CGAL::Modifier_base<HalfedgeDS> {
public:
	SurfaceCreator(const std::vector<Vector>& ptsIn,
			const std::vector<Index>& polysIn,
			const std::vector<int>& labelsIn) :
				points(ptsIn), triangles(polysIn), labels(labelsIn) {
		}
	void operator()(HalfedgeDS& hds);

private:
	const std::vector<Vector>& points;
	const std::vector<Index>& triangles;
	const std::vector<int>& labels;

};

void SurfaceCreator::operator()(HalfedgeDS& hds) {
	CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(hds, true);
	B.begin_surface(points.size(), triangles.size(), 0);

	typedef typename HalfedgeDS::Vertex Vertex;
	typedef typename Vertex::Point Point;

	for (auto pt : points)
		B.add_vertex(Point(pt.x, pt.y, pt.z));

	Face_handle face;

	for (size_t i = 0; i < triangles.size(); ++i) {
		auto& tri = triangles[i];
		if (tri[0] != tri[1] & tri[0] != tri[2] & tri[1] != tri[2]) {
			//VTK polygons can contain lines where two vertexes are identical. Forget these
			auto face = B.begin_facet();
			B.add_vertex_to_facet(tri[0]);
			B.add_vertex_to_facet(tri[1]);
			B.add_vertex_to_facet(tri[2]);
			B.end_facet();
			// The face id is size_t i.e. unsigned so we shift this to
			// positive. 1 is wall. 2,3 ... are the inlets and outlets.
			face->id() = labels[i] + 2;
		} else {
			std::cout << "Eliminated degenerate vertex: " << tri << std::endl;
		}
	}

	B.end_surface();
	B.remove_unconnected_vertices();

}

std::shared_ptr<Polyhedron> MkCGALMesh(const std::vector<Vector>& points,
		const std::vector<Index>& triangles,
		const std::vector<int>& labels) {
	SurfaceCreator modifier(points, triangles, labels);
	auto mesh = std::make_shared<Polyhedron>();
	mesh->delegate(modifier);
	return mesh;
}
