include $(MK)/header.mk

SRCS := Collision.cc \
        InletOutletCollision.cc \
        InletOutletWallCollision.cc \
        MidFluidCollision.cc \
        WallCollision.cc \
        ImplNonZeroVelocityBoundaryDensity.cc \
        ImplSimpleCollideAndStream.cc \
        ImplZeroVelocityBoundaryDensity.cc \
        ImplZeroVelocityEquilibrium.cc \
        ImplFInterpolation.cc \
        ImplGuoZhengShi.cc \
        ImplRegularised.cc \
        ImplSimpleBounceBack.cc \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
