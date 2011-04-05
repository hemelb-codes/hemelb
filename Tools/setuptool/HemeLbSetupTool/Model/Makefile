CXX = gcc

CXXFLAGS = -fno-strict-aliasing -fno-common -dynamic -arch i386 -DNDEBUG -g -O3 -arch i386 -Wno-deprecated \
	-I/Library/Frameworks/Python.framework/Versions/7.0/include \
	-I/Library/Frameworks/Python.framework/Versions/7.0/lib/python2.7/site-packages/numpy/core/include \
	-I/Library/Frameworks/Python.framework/Versions/7.0/include/python2.7 \
	-I/Library/Frameworks/Python.framework/Versions/Current/include/vtk-5.6

LDFLAGS = -L/Library/Frameworks/Python.framework/Versions/7.0/lib \
	-L/Library/Frameworks/Python.framework/Versions/7.0/lib \
	-bundle -undefined dynamic_lookup -arch i386 -g

targets = OutputGenerationHelpers.so HitList.so vtkHelp.so
all : $(targets)
$(echo )
%.cpp : %.pyx
	cython --cplus $<

.PRECIOUS : %.cpp %.o

%.so : %.o
	$(CXX) $(LDFLAGS) $^ -o $@

OutputGenerationHelpers.cpp : OutputGenerationHelpers.pyx vtkHelp.pxd HitList.pxd Flags.pxi
HitList.cpp : vtkHelp.pxd
vtkHelp.cpp : vtkHelp.pxd

.PHONY : clean
clean :
	-rm $(targets) $(addsuffix .cpp, $(basename $(targets))) $(addsuffix .o, $(basename $(targets)))

