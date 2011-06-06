#include "ConfigGenerator.h"
#include "ConfigWriter.h"
#include "Domain.h"
#include "Block.h"

// From the HemeLB tree
#include "config.h"

#include "vtkPolyDataAlgorithm.h"
#include "vtkOBBTree.h"
#include "vtkPolyData.h"

//ConfigGenerator::ConfigGenerator() {
//}

void ConfigGenerator::Execute() {
	ConfigGenerator& self = *this;

	Domain domain(self.VoxelSize, self.SurfaceSource->GetOutput()->GetBounds());

	ConfigWriter writer(self.OutputConfigFile, self.StressType,
			domain.GetBlockSize(), domain.GetBlockCounts(),
			domain.GetVoxelSize(), domain.GetOrigin());

	for (Domain::iterator block = domain.begin(); block != domain.end(); ++block) {
		// Open the BlockStarted context of the writer; this will
		// deal with flushing the state to the file (or not, in the
		//case where there are no fluid sites).
		BlockWriter* blockWriterPtr = writer.StartNextBlock();
		// Get a reference to make the syntax nicer below!
		BlockWriter& blockWriter = *blockWriterPtr;

		for (Block::iterator site = (*block).begin(); site != block->end(); ++site) {
			this->ClassifySite(*site);
			// cache the type cos it's probably slow to compute
			unsigned int type = site->GetType();
			unsigned int cfg = site->GetConfig();
			blockWriter << cfg;

			if (type == hemelb::SOLID_TYPE)
				//Solid sites, we don't do anything
				continue;

			// Increase count of fluid sites
			blockWriter.IncrementFluidSitesCount();

			if (cfg == hemelb::FLUID_TYPE)
				// Pure fluid sites don't need any more data
				continue;

			// It must be INLET/OUTLET/EDGE

			if (type == hemelb::INLET_TYPE || type == hemelb::OUTLET_TYPE) {
				for (Vector::iterator bn = site->BoundaryNormal.begin(); bn
						!= site->BoundaryNormal.end(); ++bn)
					blockWriter << *bn;
				blockWriter << site->BoundaryDistance;
			}

			if (site->IsEdge) {
				for (Vector::iterator wn = site->WallNormal.begin(); wn
						!= site->WallNormal.end(); ++wn)
					blockWriter << *wn;
				blockWriter << site->WallDistance;
			}

			for (std::vector<double>::iterator it = site->CutDistances.begin(); it
					!= site->CutDistances.end(); ++it)
				blockWriter << *it;

		}

	}
}

void ConfigGenerator::ClassifySite(Site& site) {

}

