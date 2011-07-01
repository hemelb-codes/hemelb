include $(MK)/header.mk

SRCS := CollisionVisitor.cc \
	Collision.cc \
        InletOutletCollision.cc \
        InletOutletWallCollision.cc \
        MidFluidCollision.cc \
        WallCollision.cc \
	HFunction.cc \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
