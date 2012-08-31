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

CylinderGenerator::CylinderGenerator() {
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
			double sinTheta;
			double n_i = n[i] > 0. ? n[i] : -n[i];
			double n_i2 = n_i * n_i;
			if (n_i2 > 1.) {
				sinTheta = 0.;
			} else {
				sinTheta = std::sqrt(1. - n_i2);
			}
			bounds[2 * i] = c[i] - 0.5 * h * n_i - r * sinTheta;
			bounds[2 * i + 1] = c[i] + 0.5 * h * n_i + r * sinTheta;
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

bool IsInsideCylinder(CylinderData* cyl, Vector& pt) {
	Vector x = pt - cyl->Centre;
	double z = Vector::Dot(x, cyl->Axis);
	if (z > -0.5 * cyl->Length && z < 0.5 * cyl->Length) {
		double rsq = Vector::Dot(x, x) - z * z;
		if (rsq < cyl->Radius * cyl->Radius) {
			return true;
		} else {
			return false;
		}
	} else {
		return false;
	}
}

#define TOL 1e-6

struct Hit {
	double t;
	Vector pt;
	int cellId;
	bool operator<(const Hit& rhs) const {
		// Note this goes the wrong way! It's meant to.
		return (t > rhs.t);
	}
};

#include <queue>
typedef std::priority_queue<Hit> HitList;
/**
 * Compute the intersection between the line segment and the capped cylinder.
 *
 * Site from is ALWAYS fluid
 * Site to is ALWAYS solid
 *
 * Therefore we expect exactly one intersection. Due to floating point errors,
 * we can't guarantee this. So, we search for intersections with some
 * tolerance (macro TOL) and pick the closest to fluid site (as given by the
 * parametic line coordinate t. We use a priority_queue to keep the hits in
 * order.
 */
Hit ComputeIntersection(CylinderData* cyl, std::vector<Iolet*>& iolets,
		Site& from, Site& to) {
	HitList hitList = HitList();
	/*
	 * The equation of the line segment is:
	 * x(t) = a + t (b - a)		t E (0, 1)
	 */
	Vector& n = cyl->Axis;
	double& r = cyl->Radius;
	Vector& c = cyl->Centre;
	double& h = cyl->Length;

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

		double A = (b_a.GetMagnitudeSquared() - b_aDOTn * b_aDOTn);
		double B = 2. * (Vector::Dot(a, b_a) - aDOTn * b_aDOTn);
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
					if (xDOTn >= -0.5 * h - TOL && xDOTn <= 0.5 * h + TOL) {
						// Real cylinder hit!
						Hit hit;
						hit.t = t;
						hit.pt = from.Position + b_a * t;
						hit.cellId = -1;
						hitList.push(hit);
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
	for (std::vector<Iolet*>::iterator iIt = iolets.begin(); iIt
			!= iolets.end(); ++iIt) {
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
			if (radial.GetMagnitudeSquared() < (r + TOL) * (r + TOL)) {
				// Within the cap
				Hit hit;
				hit.t = t;
				hit.pt = x;
				hit.cellId = iIt - iolets.begin();
				hitList.push(hit);
			}
		}
	}
	// We know there SHOULD be exactly one intersection. Take the closest.
	if (hitList.size() == 0)
		throw InconsistentFluidnessError(from, to, hitList.size());
	return hitList.top();
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
			neigh.IsFluidKnown = true;
			neigh.IsFluid = IsInsideCylinder(this->Cylinder, neigh.Position);
			if (neigh.IsFluid)
				neigh.CreateLinksVector();
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

			if (site.IsFluid) {
				fluid = &site;
				solid = &neigh;
				iSolid = iNeigh;
			} else {
				fluid = &neigh;
				solid = &site;
				iSolid = Neighbours::inverses[iNeigh];
			}

			// This uses special knowledge about the fact that we're working
			// on a cylinder. Since the first point is always fluid, and the
			// second always solid, we must get exactly one hit.
			Hit hit = ComputeIntersection(this->Cylinder, this->Iolets, *fluid,
					*solid);

			LinkData& link = fluid->Links[iSolid];

			// This is set in any solid case
			link.Distance = (hit.pt - fluid->Position).GetMagnitude();
			// The distance is in voxels but must be output as a fraction of
			// the lattice vector. Scale it.
			link.Distance /= Neighbours::norms[iSolid];

			if (hit.cellId < 0) {
				// -1 => we hit a wall
				link.Type = geometry::CUT_WALL;
			} else {
				// We hit an inlet or outlet
				Iolet* iolet = this->Iolets[hit.cellId];
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
