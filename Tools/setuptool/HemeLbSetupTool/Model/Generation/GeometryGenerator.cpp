//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "GeometryGenerator.h"
#include "GeometryWriter.h"

#include "Neighbours.h"
#include "Site.h"
#include "Block.h"
#include "Domain.h"

#include "Debug.h"

#include "io/formats/geometry.h"
#include <cassert>

using namespace hemelb::io::formats;

GeometryGenerator::GeometryGenerator() {
	Neighbours::Init();
}

GeometryGenerator::~GeometryGenerator() {
}

void GeometryGenerator::PreExecute() {

}

void GeometryGenerator::Execute(bool skipNonIntersectingBlocks)
		throw (GenerationError) {

	if (skipNonIntersectingBlocks) {
		throw GenerationErrorMessage(
				"Skip non intersecting blocks functionality currently not available. See ticket #651");
	}

	this->PreExecute();
	double bounds[6];
	this->ComputeBounds(bounds);
	Domain domain(this->OriginWorking, this->SiteCounts);

	GeometryWriter writer(this->OutputGeometryFile, domain.GetBlockSize(),
			domain.GetBlockCounts());

	for (BlockIterator blockIt = domain.begin(); blockIt != domain.end();
			++blockIt) {
		// Open the BlockStarted context of the writer; this will
		// deal with flushing the state to the file (or not, in the
		// case where there are no fluid sites).
		BlockWriter* blockWriterPtr = writer.StartNextBlock();
		Block& block = *blockIt;

		int side = 0; // represents whether the block is inside (-1) outside (+1) or undetermined (0)

		if (skipNonIntersectingBlocks) {
			side = this->BlockInsideOrOutsideSurface(block);
		} else { // don't use the optimisation -- check every site
			side = 0;
		}

		switch (side) {
		case 1:
			// Block is entirely outside the domain.
			// We don't have to do anything.
			break;
		case 0:
			// Block has some surface within it.
			for (SiteIterator siteIt = block.begin(); siteIt != block.end();
					++siteIt) {
				Site& site = **siteIt;
				this->ClassifySite(site);
				// here we should check site
				if (site.IsFluid) {
					blockWriterPtr->IncrementFluidSitesCount();
					WriteFluidSite(*blockWriterPtr, site);
				} else {
					WriteSolidSite(*blockWriterPtr, site);
				}

			}
			break;
		case -1:
			// Block is entirely inside the domain
			for (SiteIterator siteIt = block.begin(); siteIt != block.end();
					++siteIt) {
				Site& site = **siteIt;
				site.IsFluidKnown = true;
				site.IsFluid = true;
				site.CreateLinksVector();
				for (unsigned int link_index = 0;
						link_index < site.Links.size(); ++link_index) {
					site.Links[link_index].Type = geometry::CUT_NONE;
				}
				blockWriterPtr->IncrementFluidSitesCount();
				WriteFluidSite(*blockWriterPtr, site);
			}
			break;
		default:
			break;
		}
		blockWriterPtr->Finish();
		blockWriterPtr->Write(writer);
		delete blockWriterPtr;
	}
	writer.Close();
}

void GeometryGenerator::WriteSolidSite(BlockWriter& blockWriter, Site& site) {
	blockWriter << static_cast<unsigned int>(geometry::SOLID);
	// That's all in this case.
}

void GeometryGenerator::WriteFluidSite(BlockWriter& blockWriter, Site& site) {
	blockWriter << static_cast<unsigned int>(geometry::FLUID);

	// Iterate over the displacements of the neighbourhood
	for (unsigned int i = 0; i < Neighbours::n; ++i) {
		unsigned int cutType = site.Links[i].Type;

		if (cutType == geometry::CUT_NONE) {
			blockWriter << static_cast<unsigned int>(geometry::CUT_NONE);
		} else if (cutType == geometry::CUT_WALL
				|| cutType == geometry::CUT_INLET
				|| cutType == geometry::CUT_OUTLET) {
			blockWriter << static_cast<unsigned int>(cutType);
			if (cutType == geometry::CUT_INLET
					|| cutType == geometry::CUT_OUTLET) {
				blockWriter << static_cast<unsigned int>(site.Links[i].IoletId);
			}
			blockWriter << static_cast<float>(site.Links[i].Distance);
		} else {
			// TODO: throw some exception
			std::cout << "Unknown cut type "
					<< static_cast<unsigned int>(cutType) << " for site "
					<< site.GetIndex() << std::endl;
		}
	}

	// Indicate whether the current fluid site has a wall normal available and
	// if so write it.
	if (site.WallNormalAvailable) {
		blockWriter
				<< static_cast<unsigned int>(geometry::WALL_NORMAL_AVAILABLE);
		blockWriter << static_cast<float>(site.WallNormal[0]);
		blockWriter << static_cast<float>(site.WallNormal[1]);
		blockWriter << static_cast<float>(site.WallNormal[2]);
	} else {
		blockWriter
				<< static_cast<unsigned int>(geometry::WALL_NORMAL_NOT_AVAILABLE);
	}
}

void GeometryGenerator::ComputeAveragedNormal(Site& site) const {
	site.WallNormalAvailable = false;

	if (site.IsFluid) {
		site.WallNormal = 0.0;

		// Compute a weighted sum of the wall normals available and normalise it.
		for (unsigned neighId = 0; neighId < Neighbours::n; ++neighId) {
			LinkData& link = site.Links[neighId];
			if (link.Type == geometry::CUT_WALL) {

				assert(link.DistanceInVoxels != 0);
				double weight = 1 / link.DistanceInVoxels;
				site.WallNormal[0] += weight * link.WallNormalAtWallCut[0];
				site.WallNormal[1] += weight * link.WallNormalAtWallCut[1];
				site.WallNormal[2] += weight * link.WallNormalAtWallCut[2];
				site.WallNormalAvailable = true;
			}
		}

		// Avoid dividing by 0
		if (site.WallNormalAvailable) {
			site.WallNormal.Normalise();
		}
	}
}
