#include "ConfigGenerator.h"
#include "ConfigWriter.h"

#include "Neighbours.h"
#include "Site.h"
#include "Block.h"
#include "Domain.h"

#include "Debug.h"


#include "geometry/SiteData.h"


#include "vtkPolyDataAlgorithm.h"
#include "vtkOBBTree.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"

ConfigGenerator::ConfigGenerator() :
		ClippedSurface(NULL), IsFirstSite(true) {
	Neighbours::Init();
	this->Locator = vtkOBBTree::New();
	this->Locator->SetTolerance(1e-9);
	this->hitPoints = vtkPoints::New();
	this->hitCellIds = vtkIdList::New();
}

ConfigGenerator::~ConfigGenerator() {
	this->Locator->Delete();
	this->hitPoints->Delete();
	this->hitCellIds->Delete();
}

void ConfigGenerator::Execute() {
	// Build our locator.
	this->Locator->SetDataSet(this->ClippedSurface);
	this->Locator->BuildLocator();

	Domain domain(this->VoxelSize, this->ClippedSurface->GetBounds());

	ConfigWriter writer(this->OutputConfigFile, domain.GetBlockSize(),
			domain.GetBlockCounts(), domain.GetVoxelSize(), domain.GetOrigin());

	for (BlockIterator blockIt = domain.begin(); blockIt != domain.end();
			++blockIt) {
		// Open the BlockStarted context of the writer; this will
		// deal with flushing the state to the file (or not, in the
		// case where there are no fluid sites).
		BlockWriter* blockWriterPtr = writer.StartNextBlock();
		// Get references to make the syntax nicer below!
		BlockWriter& blockWriter = *blockWriterPtr;
		Block& block = *blockIt;

		for (SiteIterator siteIt = block.begin(); siteIt != block.end();
				++siteIt) {
			Site& site = **siteIt;
			this->ClassifySite(site);
			// cache the type cos it's probably slow to compute
			unsigned int type = site.GetType();
			uint64_t cfg = site.GetConfig();
			blockWriter << cfg;

			if (type == hemelb::geometry::SOLID_TYPE) {
				//Solid sites, we don't do anything
				continue;
			}
			// Increase count of fluid sites
			blockWriter.IncrementFluidSitesCount();

			if (cfg == hemelb::geometry::FLUID_TYPE) {
				// Pure fluid sites don't need any more data
				continue;
			}
			// It must be INLET/OUTLET/EDGE

			if (type == hemelb::geometry::INLET_TYPE || type == hemelb::geometry::OUTLET_TYPE) {
				for (Vector::iterator bn = site.BoundaryNormal.begin();
						bn != site.BoundaryNormal.end(); ++bn)
					blockWriter << *bn;
				blockWriter << site.BoundaryDistance;
			}

			if (site.IsEdge) {
				for (Vector::iterator wn = site.WallNormal.begin();
						wn != site.WallNormal.end(); ++wn)
					blockWriter << *wn;
				blockWriter << site.WallDistance;
			}
			for (std::vector<double>::iterator it = site.CutDistances.begin();
					it != site.CutDistances.end(); ++it)
				blockWriter << *it;

		}
		blockWriter.Finish();
		delete blockWriterPtr;
	}
	writer.Close();
}

bool ConfigGenerator::GetIsFluid(Site& site) {
	if (!site.IsFluidKnown) {
		if (this->Locator->InsideOrOutside(&site.Position[0]) < 0) {
			// -1 => inside surface
			site.IsFluid = true;
		} else {
			site.IsFluid = false;
		}
		site.IsFluidKnown = true;
	}
	return site.IsFluid;
}

/*
 * Perform classification of the supplied sites. Note that
 * this will alter the connected sites that have yet to be
 * classified, as we wish to examine each link only once.
 *
 * Each site must have its IsFluid and IsEdge flags set, along
 * with the CutDistances array (at appropriate indices),
 * WallDistance/Normal and BoundaryDistance/Normal.
 */
void ConfigGenerator::ClassifySite(Site& site) {

	if (!this->GetIsFluid(site)) {
		// Nothing to do for solid sites
		return;
	}

	for (NeighbourIterator neighIt = site.beginall(); neighIt != site.endall();
			++neighIt) {
		Site& neigh = *neighIt;
		if (!this->GetIsFluid(neigh)) {
			// Link to a solid site
			site.IsEdge = true;
			unsigned int iNeigh = neighIt.GetNeighbourIndex();

			this->Locator->IntersectWithLine(&site.Position[0],
					&neigh.Position[0], this->hitPoints, this->hitCellIds);
			Vector hitPoint;
			/*
			 * First hit, assign the distance (in STL units for
			 * now) to first intersection (i.e. CutDistance)
			 * and the ID of the vtkPolygon which we
			 * intersected (CutCellId).
			 */
			this->hitPoints->GetPoint(0, &hitPoint[0]);
			site.CutDistances[iNeigh] = (hitPoint - site.Position).Magnitude<
					double>();
			site.CutCellIds[iNeigh] = this->hitCellIds->GetId(0);

		}
	}

	if (!site.IsEdge)
		// Nothing more to do for simple fluid sites
		return;

	/*
	 * The CutDistances need to be fractions of the corresponding
	 * lattice vector, rather than distances in lattice units or
	 * physical units. HOWEVER, we need to know the distances to
	 * the walls and Iolets in plain lattice units below, so only
	 * scale by the VoxelSize for now and scale to vector fractions
	 * once we're done.
	 */
	for (std::vector<double>::iterator it = site.CutDistances.begin();
			it != site.CutDistances.end(); ++it)
		*it /= this->VoxelSize;

	/*
	 * Get the normals and scalars associated with the surface
	 * polygons. The scalars hold the index of the Iolet which they
	 * represent, -1 meaning they aren't an Iolet.
	 */
	vtkCellData* celldata = this->ClippedSurface->GetCellData();
	vtkDataArray* normals = celldata->GetNormals();
	vtkIntArray* ioletIds = static_cast<vtkIntArray*>(celldata->GetScalars());

	// Set this to false for now, if it's adjacent to a wall will reset later
	site.IsEdge = false;

	int hitCellId, ioletId;

	for (unsigned int iNeigh = 0; iNeigh < Neighbours::n; ++iNeigh) {
		hitCellId = site.CutCellIds[iNeigh];
		/* The site.CutCellsIds array is initialised to -1 and then
		 * updated with the index of the first vtkPolygon they
		 * intersect.
		 */
		if (hitCellId == -1) {
			// We didn't hit in this direction.
			continue;
		}

		ioletId = ioletIds->GetValue(hitCellId);
		if (ioletId >= 0) {
			// It's an iolet
			if (site.CutDistances[iNeigh] < site.BoundaryDistance) {
				// It is the closest yet
				site.AdjacentIolet = this->Iolets[ioletId];
				// TODO: confirm this should be in lattice units
				site.BoundaryDistance = site.CutDistances[iNeigh];
				site.BoundaryNormal = site.AdjacentIolet->Normal;
				site.BoundaryId = site.AdjacentIolet->Id;
			}

		} else {
			// It's wall
			site.IsEdge = true;
			if (site.CutDistances[iNeigh] < site.WallDistance) {
				// If it's the closest point yet, store it
				site.WallDistance = site.CutDistances[iNeigh];
				normals->GetTuple(hitCellId, &site.WallNormal[0]);
			}
		}

	}
	/*
	 * Scale to be fractions of lattice vectors instead of
	 * distances in lattice units; see above for more.
	 */
	for (unsigned int i = 0; i < Neighbours::n; ++i)
		site.CutDistances[i] /= Neighbours::norms[i];
}

