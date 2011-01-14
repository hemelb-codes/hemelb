include $(MK)/header.mk

TARGETS := libParmetis.$(LIBEXT)

SRCS := comm.c util.c debug.c setup.c grsetup.c timer.c \
       node_refine.c initmsection.c order.c \
       xyzpart.c pspases.c frename.c \
       iintsort.c iidxsort.c ikeysort.c ikeyvalsort.c \
       kmetis.c gkmetis.c ometis.c \
       initpart.c match.c \
       kwayfm.c kwayrefine.c kwaybalance.c \
       remap.c stat.c fpqueue.c \
       ametis.c rmetis.c lmatch.c initbalance.c \
       mdiffusion.c diffutil.c wave.c \
       csrmatch.c redomylink.c balancemylink.c \
       selectq.c akwayfm.c serial.c move.c \
       mmetis.c mesh.c memory.c weird.c backcompat.c 

$(TARGETS)_DEPS = $(subst .c,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
