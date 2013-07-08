//#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/double.h>
#include <CGAL/Random.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

//typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::Point_3 PointCGAL;
typedef Kernel::Plane_3 PlaneCGAL;
typedef Kernel::Vector_3 VectorCGAL;
typedef Kernel::Segment_3 SegmentCGAL;
typedef Kernel::Triangle_3 TriangleCGAL;
typedef Kernel::Ray_3 RayCGAL;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef Tree::Object_and_primitive_id Object_and_primitive_id;
typedef Polyhedron::Face_handle FacehandleCGAL;
typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef std::pair<Object_and_primitive_id, double> Object_Primitive_and_distance;
