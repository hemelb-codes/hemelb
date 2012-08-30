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
	$(CC) $(CFLAGS)  $< -o $@ -lglut -lGLU -lGL 


.PHONY : clean

clean :
	-rm $(target) $(objects)
