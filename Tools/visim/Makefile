target = visim

CC = g++
CFLAGS = -O3

objects = visualize_images.o
LIBS = -framework OpenGL -framework GLUT

$(target) : $(objects)
	$(CC) $(LIBS) $(objects) -o $@


.PHONY : clean

clean :
	-rm $(target) $(objects)
