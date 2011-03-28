CXX = gcc

CXXFLAGS = -fno-strict-aliasing -fno-common -dynamic -arch i386 -DNDEBUG -g -O3 -arch i386 \
	-I/Library/Frameworks/Python.framework/Versions/7.0/include \
	-I/Library/Frameworks/Python.framework/Versions/7.0/lib/python2.7/site-packages/numpy/core/include \
	-I/Library/Frameworks/Python.framework/Versions/7.0/include/python2.7 \
	-I/Library/Frameworks/Python.framework/Versions/Current/include/vtk-5.6

LDFLAGS = -L/Library/Frameworks/Python.framework/Versions/7.0/lib \
	-L/Library/Frameworks/Python.framework/Versions/7.0/lib \
	-bundle -undefined dynamic_lookup -arch i386 -g

targets = OutputGenerationHelpers.so HitList.so
all : $(targets)
$(echo )
%.cpp : %.pyx
	cython --cplus $<

.PRECIOUS : %.cpp %.o

%.so : %.o
	$(CXX) $(LDFLAGS) $^ -o $@

OutputGenerationHelpers.cpp : OutputGenerationHelpers.pyx Flags.pxi
# OutputGenerationHelpers.o : OutputGenerationHelpers.c
# 	gcc   -c OutputGenerationHelpers.c -o build/temp.macosx-10.5-i386-2.7/OutputGenerationHelpers.o

# OutputGenerationHelpers.so : OutputGenerationHelpers.o
# 	$(CXX) $(LDFLAGS) $^ -o $@

# HitList.so : HitList.o
# 	$(CXX) $(LDFLAGS) $^ -o $@

.PHONY : clean
clean :
	-rm $(targets) $(addsuffix .cpp, $(basename $(targets))) $(addsuffix .o, $(basename $(targets)))

