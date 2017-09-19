#ifndef HEMELBSETUPTOOL_CGALSEARCH_H
#define HEMELBSETUPTOOL_CGALSEARCH_H

#include "Cgal.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
//#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
namespace {
	//typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel, CgalPolyhedron> Primitive;
	typedef CGAL::AABB_face_graph_triangle_primitive<CgalPolyhedron> Primitive;
	typedef CGAL::AABB_traits<Kernel, Primitive> AABBTraits;
}

typedef CGAL::AABB_tree<AABBTraits> CgalSearchTree;

#endif
