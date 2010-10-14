include $(MK)/header.mk

SRCS := Collision.cc \
        ImplFInterpolation.cc \
        ImplGuoZhengShi.cc \
        ImplNonZeroVelocityBoundaryDensity.cc \
        ImplRegularised.cc \
        ImplSimpleBounceBack.cc \
        ImplSimpleCollideAndStream.cc \
        ImplZeroVelocityBoundaryDensity.cc \
        ImplZeroVelocityEquilibrium.cc 

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
