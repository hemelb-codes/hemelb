#include "GeometryGenerator.h"
#include "GeometryWriter.h"

#include "Neighbours.h"
#include "Site.h"
#include "Block.h"
#include "Domain.h"

#include "Debug.h"

#include "io/formats/geometry.h"

using namespace hemelb::io::formats;

GeometryGenerator::GeometryGenerator() {
	Neighbours::Init();
}

GeometryGenerator::~GeometryGenerator() {
}

void GeometryGenerator::PreExecute() {

}

void GeometryGenerator::Execute() throw (GenerationError) {
	this->PreExecute();
	double bounds[6];
	this->ComputeBounds(bounds);
	Domain domain(this->VoxelSizeMetres, bounds);

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
