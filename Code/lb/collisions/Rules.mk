include $(MK)/header.mk

SUBDIRS := implementations

TARGETS := Collision.o \
	   CollisionVisitor.o \
	   InletOutletCollision.o \
	   InletOutletWallCollision.o \
	   MidFluidCollision.o \
	   WallCollision.o \
	   HFunction.o \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
