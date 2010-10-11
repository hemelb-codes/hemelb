include $(MK)/header.mk

SRCS := Collision.cc ImplNonZeroVelocityBoundaryDensity.cc ImplSimpleCollideAndStream.cc \
  ImplZeroVelocityBoundaryDensity.cc ImplZeroVelocityEquilibrium.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
