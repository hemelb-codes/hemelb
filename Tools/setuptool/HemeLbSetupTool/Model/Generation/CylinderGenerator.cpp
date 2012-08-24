#include "CylinderGenerator.h"
#include "GeometryWriter.h"

#include "Neighbours.h"
#include "Site.h"
#include "Block.h"
#include "Domain.h"
#include "InconsistentFluidnessError.h"

#include "Debug.h"

#include "io/formats/geometry.h"

#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"

#include <sstream>
#include <cmath>

using namespace hemelb::io::formats;

CylinderGenerator::CylinderGenerator() :
	ClippedSurface(NULL) {
	Neighbours::Init();
	this->hitPoints = vtkPoints::New();
	this->hitCellIds = vtkIdList::New();
	this->Cylinder = new CylinderData;
}

CylinderGenerator::~CylinderGenerator() {
	this->hitPoints->Delete();
	this->hitCellIds->Delete();
	delete this->Cylinder;
}

void CylinderGenerator::Execute() throw (GenerationError) {
	/*
	 * Compute the approximate bounds of the cylinder
	 */
	double* bounds = new double[6];
	{
		Vector& c = this->Cylinder->Centre;
		Vector& n = this->Cylinder->Axis;
		double& r = this->Cylinder->Radius;
		double& h = this->Cylinder->Length;

		for (int i = 0; i < 3; ++i) {
			double n_i = n[i] > 0. ? n[i] : -n[i];
			bounds[2 * i] = c[i] - 0.5 * h * n_i - r;
			bounds[2 * i + 1] = c[i] + 0.5 * h * n_i + r;
		}
	}
	Domain domain(this->VoxelSizeMetres, bounds);
	delete[] bounds;

	GeometryWriter writer(this->OutputGeometryFile, domain.GetBlockSize(),
			domain.GetBlockCounts(), domain.GetVoxelSizeMetres(),
			domain.GetOriginMetres());

	for (BlockIterator blockIt = domain.begin(); blockIt != domain.end(); ++blockIt) {
		// Open the BlockStarted context of the writer; this will
		// deal with flushing the state to the file (or not, in the
		// case where there are no fluid sites).
		BlockWriter* blockWriterPtr = writer.StartNextBlock();
		Block& block = *blockIt;

		for (SiteIterator siteIt = block.begin(); siteIt != block.end(); ++siteIt) {
			Site& site = **siteIt;
			/*
			 * ClassifySite expects to be given a site of known fluidness.
			 * The constructor of Block will ensure that all sites at the edge
			 * of the Domain will be set to solid. The iterators here ensure
			 * that we start with the site at (0,0,0).
			 */
			this->ClassifySite(site);

			if (site.IsFluid) {
				blockWriterPtr->IncrementFluidSitesCount();
				WriteFluidSite(*blockWriterPtr, site);
			} else {
				WriteSolidSite(*blockWriterPtr, site);
			}
		}
		blockWriterPtr->Finish();
		blockWriterPtr->Write(writer);
		delete blockWriterPtr;
	}
	writer.Close();
}

void CylinderGenerator::WriteSolidSite(BlockWriter& blockWriter, Site& site) {
	blockWriter << static_cast<unsigned int> (geometry::SOLID);
	// That's all in this case.
}

void CylinderGenerator::WriteFluidSite(BlockWriter& blockWriter, Site& site) {
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
void CylinderGenerator::ClassifySite(Site& site) {

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

			// The value associated with that cell, which identifies what was hit.
			int ioletId = this->hitCellIds->GetId(iHit);

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

int CylinderGenerator::ComputeIntersections(Site& from, Site& to) {
	int nHits = 0;
	this->hitPoints->SetNumberOfPoints(0);
	this->hitCellIds->SetNumberOfIds(0);

	/*
	 * The equation of the line segment is:
	 * x(t) = a + t (b - a)		t E (0, 1)
	 */
	Vector& n = this->Cylinder->Axis;
	double& r = this->Cylinder->Radius;
	Vector& c = this->Cylinder->Centre;
	double& h = this->Cylinder->Length;

	{
		/*
		 * The equation determining intersection with an INFINITE cylinder of
		 * radius r is:
		 * [(b-a)^2 - ((b-a).n)^2] t^2 + [2 a.(b - a) - 2 (a.n)((b-a).n)] t + [a^2 - (a.n)^2 - r^2] = 0
		 * So first compute the coefficients and then the discriminant of the eqn
		 */
		Vector a = from.Position - c;
		Vector b = to.Position - c;
		Vector b_a = b - a;

		double b_aDOTn = Vector::Dot(b_a, n);
		double aDOTn = Vector::Dot(a, n);

		double A = (b_a.GetMagnitudeSquared() - b_aDOTn);
		double B = 2 * Vector::Dot(a, b_a) - 2 * aDOTn * b_aDOTn;
		double C = a.GetMagnitudeSquared() - aDOTn * aDOTn - r * r;

		double discriminant = B * B - 4 * A * C;

		if (discriminant < 0.) {
			// No real solutions.
		} else if (discriminant == 0) {
			// Exactly one solution, i.e. line segment just brushes the cylinder.
			// This means the line must be outside the cylinder everywhere else,
			// so we will count this as no intersection.
		} else {
			// Two real solutions. So we have two intersections between the
			// infinite line and infinite cylinder.
			// If t outside (0,1), then the intersection isn't on the line segment
			// we care about.
			std::vector<double> ts(2);
			ts[0] = (-B + std::sqrt(discriminant)) / (2 * A);
			ts[1] = (-B - std::sqrt(discriminant)) / (2 * A);
			for (std::vector<double>::iterator tIt = ts.begin(); tIt
					!= ts.end(); ++tIt) {
				double t = *tIt;
				if (t > 0. && t < 1.) {
					// Hit in part of line we care about.
					// Now check if it's on the finite cylinder. This requires that
					// x.n E (-h/2, h/2)
					double xDOTn = aDOTn + t * b_aDOTn;
					if (xDOTn >= -0.5 * h && xDOTn <= 0.5 * h) {
						// Real cylinder hit!
						Vector hit = a + b_a * t;
						++nHits;
						this->hitPoints->InsertNextPoint(hit.x, hit.y, hit.z);
						this->hitCellIds->InsertNextId(-1);
					}
				}
			}
		}
	}
	/*
	 * Now we want to look for intersections with the capping planes.
	 */
	Vector& a = from.Position;
	Vector& b = to.Position;
	Vector b_a = b - a;
	for (std::vector<Iolet*>::iterator iIt = this->Iolets.begin(); iIt
			!= this->Iolets.end(); ++iIt) {
		Iolet* iolet = *iIt;
		/*
		 * Plane equation is x.p = q.p (p = plane normal, q = point on plane)
		 * Line is x = a + t(b-a)
		 */
		Vector& q = iolet->Centre;
		Vector& p = iolet->Normal;

		double t = Vector::Dot(q - a, p) / Vector::Dot(b_a, p);
		if (t > 0. && t < 1.) {
			// Intersection within the line segment. Now check within cap.
			Vector x = a + b_a * t;
			Vector x_c = x - c;
			double x_cDOTn = Vector::Dot(x_c, n);
			Vector radial = x_c - n * x_cDOTn;
			if (radial.GetMagnitudeSquared() < r * r) {
				// Within the cap
				++nHits;
				this->hitPoints->InsertNextPoint(x.x, x.y, x.z);
				this->hitCellIds->InsertNextId(iIt - this->Iolets.begin());
			}
		}
	}
}
