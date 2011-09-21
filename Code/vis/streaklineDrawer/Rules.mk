include $(MK)/header.mk

SRCS :=	Particle.cc \
	Particles.cc \
	StreaklineDrawer.cc \
	VelocityField.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk