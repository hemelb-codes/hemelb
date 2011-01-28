display:
	g++ displayClient.cc -O3 -lGLU -lglut -o displayClient

mac-display:
	g++ -O3 -framework OpenGL -framework GLUT -framework Foundation -o displayClient displayClient.cc eVizRLEUtil.c -DMAC_GLUT

display-debug:
	g++ displayClient.cc -g -lGLU -lglut -o displayClient

clean:
	-rm displayClient
