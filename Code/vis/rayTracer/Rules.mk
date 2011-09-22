include $(MK)/header.mk

SRCS :=	ClusterNormal.cc \
	ClusterWithWallNormals.cc \
	HSLToRGBConverter.cc \
	RayDataEnhanced.cc \
	RayDataNormal.cc \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
