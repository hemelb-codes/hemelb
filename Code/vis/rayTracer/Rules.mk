include $(MK)/header.mk

SRCS :=	BlockTraverser.cc \
	BlockTraverserWithVisitedBlockTracker.cc \
	Cluster.cc \
	ClusterWithWallNormals.cc \
	ClusterBuilder.cc \
	ClusterBuilderWithWallNormals.cc \
	ClusterRayTracer.cc \
	ClusterRayTracerEnhanced.cc \
	ClusterTraverser.cc \
	Ray.cc \
	RayEnhanced.cc \
	RayTracer.cc \
	RayTracerEnhanced.cc \
	VolumeTraverser.cc \
	SiteTraverser.cc \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
