target = visim

CC = g++
CFLAGS = -O4 -L/usr/include/GT/ -lglut -lGLU

objects = visualize_images.o

$(target) : $(objects)
	$(CC) $(CFLAGS) $(objects) -o $@


.PHONY : clean

clean :
	-rm $(target) $(objects)
