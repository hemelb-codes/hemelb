include $(MK)/header.mk

SRCS :=	BlockTraverser.cc \
	Cluster.cc \
	ClusterWithWallNormals.cc \
	ClusterBuilder.cc \
	ClusterBuilderWithWallNormals.cc \
	RayTracer.cc \
	RayTracerWithLighting.cc \
	VolumeTraverser.cc \
	SiteTraverser.cc \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
