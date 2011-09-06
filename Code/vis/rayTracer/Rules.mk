include $(MK)/header.mk

SRCS :=	BlockTraverser.cc \
	BlockTraverserWithVisitedBlockTracker.cc \
	ClusterNormal.cc \
	ClusterWithWallNormals.cc \
	HSLToRGBConverter.cc \
	RayDataEnhanced.cc \
	RayDataNormal.cc \
	VolumeTraverser.cc \
	SiteTraverser.cc \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
