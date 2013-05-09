//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "PolyDataGenerator.h"

#include "Neighbours.h"
#include "Site.h"
#include "InconsistentFluidnessError.h"

#include "Debug.h"

#include "io/formats/geometry.h"

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

#include <iostream>

using namespace hemelb::io::formats;

PolyDataGenerator::PolyDataGenerator() :
		GeometryGenerator(), ClippedSurface(NULL) {
	this->Locator = vtkOBBTree::New();
	//this->Locator->SetNumberOfCellsPerNode(32); // the default
	this->Locator->SetTolerance(1e-9);
	this->hitPoints = vtkPoints::New();
	this->hitCellIds = vtkIdList::New();	
}

PolyDataGenerator::~PolyDataGenerator() {
	this->Locator->Delete();
	this->hitPoints->Delete();
	this->hitCellIds->Delete();
}

void PolyDataGenerator::ComputeBounds(double bounds[6]) const {
	this->ClippedSurface->GetBounds(bounds);
	//	for (int i=0;i<6;i++) {
	//  cout << bounds[i] << endl;
	//}
}

void PolyDataGenerator::PreExecute(void) {
	// Build our locator.
	this->Locator->SetDataSet(this->ClippedSurface);
	this->Locator->BuildLocator();

	/*
	 * Get the scalars associated with the surface polygons. The scalars hold
	 * the index of the Iolet which they represent, -1 meaning they aren't an
	 * Iolet.
	 */
	this->IoletIdArray = vtkIntArray::SafeDownCast(
			this->ClippedSurface->GetCellData()->GetScalars());
	if (this->IoletIdArray == NULL) {
		throw GenerationErrorMessage(
				"Error getting Iolet ID array from clipped surface");
	}
	
}

/*
 * Given a site with known fluidness, examine the links to not-yet-visited
 * neighbouring sites. If the neighbours have unknown fluidness, set that.
 *
 * Since we wish to examine each link only once, this will set link properties
 * from neighbour => site as well as site => neighbour
 *
 */
void PolyDataGenerator::ClassifySite(Site& site) {

	for (LaterNeighbourIterator neighIt = site.begin(); neighIt != site.end();
			++neighIt) {
		Site& neigh = *neighIt;
		unsigned int iNeigh = neighIt.GetNeighbourIndex();
		vtkIdType nHits;

		if (!neigh.IsFluidKnown) {
			// Neighbour unknown, must always intersect
			nHits = this->ComputeIntersections(site, neigh);

			if (nHits % 2 == 0) {
				// Even # hits, hence neigh has same type as site
				neigh.IsFluid = site.IsFluid;
			} else {
				// Odd # hits, neigh is opposite type to site
				neigh.IsFluid = !site.IsFluid;
			}

			if (neigh.IsFluid)
				neigh.CreateLinksVector();

			neigh.IsFluidKnown = true;
		} else {
			// We know the fluidness of neigh, maybe don't need to intersect
			if (site.IsFluid != neigh.IsFluid) {
				// Only in the case of difference must we intersect.
				nHits = this->ComputeIntersections(site, neigh);
				if (nHits % 2 != 1) {
					throw InconsistentFluidnessError(site, neigh, nHits);
				}
			}
		}

		// Four cases: fluid-fluid, solid-solid, fluid-solid and solid-fluid.
		// Will handle the last two together.
		if (site.IsFluid == neigh.IsFluid) {
			if (site.IsFluid) {
				// Fluid-fluid, must set CUT_NONE for both
				site.Links[iNeigh].Type = geometry::CUT_NONE;
				neigh.Links[Neighbours::inverses[iNeigh]].Type =
						geometry::CUT_NONE;
			} else {
				// solid-solid, nothing to do.
			}
		} else {
			// They differ, figure out which is fluid and which is solid.
			Site* fluid;
			Site* solid;

			// Index of the solid site from the fluid site.
			int iSolid;
			// Index of the point in this->hitPoints we're considering as the
			// hit of interest (i.e. the one closest to the fluid site).
			int iHit;

			if (site.IsFluid) {
				fluid = &site;
				solid = &neigh;
				iSolid = iNeigh;
				iHit = 0;
			} else {
				fluid = &neigh;
				solid = &site;
				iSolid = Neighbours::inverses[iNeigh];
				iHit = nHits - 1;
			}

			Vector hitPoint;
			this->hitPoints->GetPoint(iHit, &hitPoint[0]);
			LinkData& link = fluid->Links[iSolid];

			// This is set in any solid case
			float distanceInVoxels =
					(hitPoint - fluid->Position).GetMagnitude();
			// The distance is in voxels but must be output as a fraction of
			// the lattice vector. Scale it.
			link.Distance = distanceInVoxels / Neighbours::norms[iSolid];

			// The index of the cell in the vtkPolyData that was hit
			int hitCellId = this->hitCellIds->GetId(iHit);
			// The value associated with that cell, which identifies what was hit.
			int ioletId = this->IoletIdArray->GetValue(hitCellId);

			if (ioletId < 0) {
				// -1 => we hit a wall
				link.Type = geometry::CUT_WALL;
			} else {
				// We hit an inlet or outlet
				Iolet* iolet = this->Iolets[ioletId];
				if (iolet->IsInlet) {
					link.Type = geometry::CUT_INLET;
				} else {
					link.Type = geometry::CUT_OUTLET;
				}
				// Set the Id
				link.IoletId = iolet->Id;
			}

			// If this link intersected the wall, store the normal of the cell we hit and the distance to it.
			if (link.Type == geometry::CUT_WALL) {
				double* normal =
						this->Locator->GetDataSet()->GetCellData()->GetNormals()->GetTuple3(
								hitCellId);
				link.WallNormalAtWallCut = Vector(normal[0], normal[1],
						normal[2]);
				link.DistanceInVoxels = distanceInVoxels;
			}
		}
	}

	// If there's enough information available, an approximation of the wall normal will be computed for this fluid site.
	this->ComputeAveragedNormal(site);
}

int PolyDataGenerator::ComputeIntersections(Site& from, Site& to) {
	this->Locator->IntersectWithLine(&from.Position[0], &to.Position[0],
			this->hitPoints, this->hitCellIds);
	
	//cout << this->hitPoints->GetNumberOfPoints() << endl;
	std::ifstream     in;
	const char* infile = "working_cylinder_clipped.off"; 
        in.open(infile);
	std::istream* p_in = &in;
	Polyhedron P;
	(*p_in) >> P;
	if (!*p_in) {
	  std::cerr << "error: cannot open file"<< std::endl;
	  exit( 1);
	}
	this->ClippedCGALSurface = &P;
	Tree tree(this->ClippedCGALSurface->facets_begin(),this->ClippedCGALSurface->facets_end()); 
	this->AABBtree =  &tree;
	//cout << this->ClippedCGALSurface->size_of_facets() << endl;
	
	//cout << this->AABBtree->size() << endl;
	//cout << &this->AABBtree << endl;
	//cout << &this->ClippedCGALSurface << endl;
	

	// Vertex_iteratorCGAL vi = P.vertices_begin();
	// PointCGAL p = vi->point();
	// double minx = p.x();
	// double miny = p.y();
	// double minz = p.z();
	// double maxx = p.x();
	// double maxy = p.y();
	// double maxz = p.z();
	// for ( ; vi != P.vertices_end() ; ++vi) {
	//   p = vi->point();
	//   if ( p.x() < minx)
        //     minx = p.x();
	//   if ( p.y() < miny)
        //     miny = p.y();
	//   if ( p.z() < minz)
        //     minz = p.z();
	//   if ( p.x() > maxx)
        //     maxx = p.x();
	//   if ( p.y() > maxy)
        //     maxy = p.y();
	//   if ( p.z() > maxz)
        //     maxz = p.z();
	// }
	//cout << minx << " " << miny << " " << minz << " " << endl;
	//Cout << maxx << " " << maxy << " " << maxz << " " << endl;
	//PointCGAL p1(-9.5, -4.5, -120.5);
	//PointCGAL p2(-8.5, -3.5, -119.5);
	//SegmentCGAL segment_query(p1,p2);
	//cout << tree.number_of_intersected_primitives(segment_query) << endl;
	PointCGAL p(from.Position[0], from.Position[1], from.Position[2]);
	PointCGAL q(to.Position[0], to.Position[1], to.Position[2]);
	SegmentCGAL segment_query(p,q);
	//cout << from.Position[0] << " "<< from.Position[1] << " " << from.Position[2] << endl; 
	//cout << to.Position[0] << " " << to.Position[1]<< " " << to.Position[2] << endl; 
	//cout << from.Position[0] - to.Position[0] << " " << from.Position[1] - to.Position[1]<< " " << from.Position[2] - to.Position[2] << endl; 
	//cout << this->AABBtree->size() << endl;
	//cout << "and" << endl;
	//cout << this->ClippedCGALSurface->size_of_facets() << endl;
	//cout << &this->AABBtree << endl;
	//cout << "and" << endl;
	//cout << &this->ClippedCGALSurface << endl;
	int j = this->AABBtree->number_of_intersected_primitives(segment_query);
	PointCGAL point1;
	PointCGAL point2;
	if (j != this->hitPoints->GetNumberOfPoints()){
	  std::vector<Object_and_primitive_id> intersections;
	  this->AABBtree->all_intersections(segment_query,std::back_inserter(intersections));
	  if (CGAL::assign(point1, intersections[0].first) && CGAL::assign(point2, intersections[1].first)) {
	    //cout << "Point1 " << point1.x() << " " << point1.y() << " " << point1.z() << endl;
	    //cout << "Point2 " << point2.x() << " " << point2.y() << " " << point2.z() << endl;
	    CGAL::Comparison_result pointsidentical = CGAL::compare_xyz(point1,point2);
	    cout << pointsidentical << endl;
	    //cout << j << " is not equal to "<< this->hitPoints->GetNumberOfPoints() << endl;
	    //cout << "From " << from.Position[0] << " " << from.Position[1] << " " << from.Position[2] << endl;
	    //cout << "To " << to.Position[0] << " " << to.Position[1] << " " << to.Position[2] << endl;
	  }
	  else {
	    cout << "Not just points" << endl;
	  }
	  
	}
	return this->hitPoints->GetNumberOfPoints();
}

// Function to be called on intersecting leaf nodes of the two OBB trees.
// Final void pointer is a pointer to an int, namely the count of the number
// of intersections found so far, which is incremented.
int IntersectingLeafCounter(vtkOBBNode* polyNode, vtkOBBNode* cubeNode,
		vtkMatrix4x4* transform, void *ptr_to_intersection_count) {
	int &intersection_count = *static_cast<int*>(ptr_to_intersection_count);
	intersection_count++;
}

int PolyDataGenerator::BlockInsideOrOutsideSurface(const Block &block) {
	// Create an OBB tree for the block
	vtkOBBTree *blockSlightlyLargerOBBTree = block.CreateOBBTreeModel(1.0);

	// Count the number of domain OBB leaf nodes that intersect the single
	// node created for the block.
	int intersection_count = 0;
	Locator->IntersectWithOBBTree(blockSlightlyLargerOBBTree, NULL,
			IntersectingLeafCounter, static_cast<void*>(&intersection_count));
	// Delete the underlying polydata
	blockSlightlyLargerOBBTree->GetDataSet()->Delete();
	// And the OBBTree itself
	blockSlightlyLargerOBBTree->Delete();

	if (intersection_count == 0) {
		// either entirely inside or entirely outside
		double middlePosition[3];
		middlePosition[0] = block.Middle().Position[0];
		middlePosition[1] = block.Middle().Position[1];
		middlePosition[2] = block.Middle().Position[2];
		return Locator->InsideOrOutside(middlePosition);
	}
	return 0;
}
