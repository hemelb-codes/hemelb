include $(MK)/header.mk

SRCS := BlockTraverser.cc \
        Cluster.cc \
        ClusterBuilder.cc \
        RayTracer.cc \
        VolumeTraverser.cc \
        SiteTraverser.cc \
        RayPixel.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
