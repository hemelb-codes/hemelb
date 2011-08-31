include $(MK)/header.mk

SRCS :=	BlockTraverser.cc \
	BlockTraverserWithVisitedBlockTracker.cc \
	ClusterWithWallNormals.cc \
	Ray.cc \
	RayEnhanced.cc \
	VolumeTraverser.cc \
	SiteTraverser.cc \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
