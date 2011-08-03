include $(MK)/header.mk

SRCS :=	BlockIterator.cc \
	Location.cc \
	RectangularIterator.cc \
	ClusterBuilder.cc \
	RayTracer.cc \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
