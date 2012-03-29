#include "GeometryGenerator.h"
#include "GeometryWriter.h"

#include "Neighbours.h"
#include "Site.h"
#include "Block.h"
#include "Domain.h"

#include "Debug.h"

#include "io/formats/geometry.h"

#include "vtkPolyDataAlgorithm.h"
#include "vtkOBBTree.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"

#include <sstream>

using namespace hemelb::io::formats;

void FormatSite(std::ostringstream& msg, const Site& site) {
	msg << "site index " << site.GetIndex() << ", position " << site.Position
			<< ", which is ";
	if (site.IsFluid)
		msg << "fluid";
	else
		msg << "solid";
}

InconsistentFluidnessError::InconsistentFluidnessError(const Site& s1,
		const Site& s2, const int nHits) :
	s1(s1), s2(s2), nHits(nHits) {
	std::ostringstream msgStream;
	msgStream << "Inconsistent fluidness detected between ";
	FormatSite(msgStream, s1);
	msgStream << " and ";
	FormatSite(msgStream, s2);
	msgStream << " but found " << nHits << " intersections with the surface.";
	this->msg = msgStream.str();
}


const char* InconsistentFluidnessError::what() const throw () {
	return this->msg.c_str();
}

GeometryGenerator::GeometryGenerator() :
	ClippedSurface(NULL) {
	Neighbours::Init();
	this->Locator = vtkOBBTree::New();
	this->Locator->SetTolerance(1e-9);
	this->hitPoints = vtkPoints::New();
	this->hitCellIds = vtkIdList::New();
}

GeometryGenerator::~GeometryGenerator() {
	this->Locator->Delete();
	this->hitPoints->Delete();
	this->hitCellIds->Delete();
}

void GeometryGenerator::Execute() throw (GenerationError) {
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

	Domain domain(this->VoxelSizeMetres, this->ClippedSurface->GetBounds());

	GeometryWriter writer(this->OutputGeometryFile, domain.GetBlockSize(),
			domain.GetBlockCounts(), domain.GetVoxelSizeMetres(),
			domain.GetOriginMetres());

	for (BlockIterator blockIt = domain.begin(); blockIt != domain.end(); ++blockIt) {
		// Open the BlockStarted context of the writer; this will
		// deal with flushing the state to the file (or not, in the
		// case where there are no fluid sites).
		BlockWriter& blockWriter = writer.StartNextBlock();
		Block& block = *blockIt;

		for (SiteIterator siteIt = block.begin(); siteIt != block.end(); ++siteIt) {
			Site& site = *siteIt;
			/*
			 * ClassifySite expects to be given a site of known fluidness.
			 * The constructor of Block will ensure that all sites at the edge
			 * of the Domain will be set to solid. The iterators here ensure
			 * that we start with the site at (0,0,0).
			 */
			this->ClassifySite(site);

			if (site.IsFluid) {
				blockWriter.IncrementFluidSitesCount();
				WriteFluidSite(blockWriter, site);
			} else {
				WriteSolidSite(blockWriter, site);
			}
		}
		blockWriter.Finish();
		blockWriter.Write(writer);
	}
	writer.Close();
}

void GeometryGenerator::WriteSolidSite(BlockWriter& blockWriter, Site& site) {
	blockWriter << static_cast<unsigned int> (geometry::SOLID);
	// That's all in this case.
}

void GeometryGenerator::WriteFluidSite(BlockWriter& blockWriter, Site& site) {
	blockWriter << static_cast<unsigned int> (geometry::FLUID);

	// Iterate over the displacements of the neighbourhood
	for (unsigned int i = 0; i < Neighbours::n; ++i) {
		unsigned int cutType = site.Links[i].Type;

		if (cutType == geometry::CUT_NONE) {
			blockWriter << static_cast<unsigned int> (geometry::CUT_NONE);
		} else if (cutType == geometry::CUT_WALL || cutType
				== geometry::CUT_INLET || cutType == geometry::CUT_OUTLET) {
			blockWriter << static_cast<unsigned int> (cutType);
			if (cutType == geometry::CUT_INLET || cutType
					== geometry::CUT_OUTLET) {
				blockWriter
						<< static_cast<unsigned int> (site.Links[i].IoletId);
			}
			blockWriter << static_cast<float> (site.Links[i].Distance);
		} else {
			// TODO: throw some exception
			std::cout << "Unknown cut type "
					<< static_cast<unsigned int> (cutType) << " for site "
					<< site.GetIndex() << std::endl;
		}
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
void GeometryGenerator::ClassifySite(Site& site) {

	for (LaterNeighbourIterator neighIt = site.begin(); neighIt != site.end(); ++neighIt) {
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
				neigh.Links[Neighbours::inverses[iNeigh]].Type
						= geometry::CUT_NONE;
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
			link.Distance = (hitPoint - fluid->Position).GetMagnitude();
			// The distance is in voxels but must be output as a fraction of
			// the lattice vector. Scale it.
			link.Distance /= Neighbours::norms[iSolid];

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
		}

	}
}

int GeometryGenerator::ComputeIntersections(Site& from, Site& to) {
	this->Locator->IntersectWithLine(&from.Position[0], &to.Position[0],
			this->hitPoints, this->hitCellIds);
	return this->hitPoints->GetNumberOfPoints();
}
