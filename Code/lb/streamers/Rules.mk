include $(MK)/header.mk

TARGETS := Collision.o \
           InletOutletCollision.o \
           InletOutletWallCollision.o \
           MidFluidCollision.o \
           WallCollision.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
