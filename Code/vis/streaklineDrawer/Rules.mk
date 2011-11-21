include $(MK)/header.mk

SRCS := NeighbouringProcessor.cc \
        Particle.cc \
        ParticleManager.cc \
        StreaklineDrawer.cc \
        StreakPixel.cc \
        VelocityField.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk