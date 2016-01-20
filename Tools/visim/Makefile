
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
target = visim

CC = g++
ifeq ($(HEMELB_MACHINE),OSX)
CFLAGS = -framework OpenGL -framework GLUT -framework Foundation -O4 -DHEMELB_CFG_ON_OSX
else
CFLAGS = -O4 -L/usr/include/GT/ -lglut -lGLU -lGL
endif

objects = visualize_images.o

src =  visualize_images.cc
$(target) : $(src)
	$(CC) $(CFLAGS)  $< -o $@ 


.PHONY : clean

clean :
	-rm $(target) $(objects)
