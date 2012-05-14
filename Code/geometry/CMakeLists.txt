add_library(
	hemelb_geometry BlockTraverser.cc BlockTraverserWithVisitedBlockTracker.cc 
	GeometryReader.cc Decomposition.cc LatticeData.cc SiteData.cc SiteTraverser.cc VolumeTraverser.cc
	Block.cc decomposition/BasicDecomposition.cc decomposition/OptimisedDecomposition.cc
	)
