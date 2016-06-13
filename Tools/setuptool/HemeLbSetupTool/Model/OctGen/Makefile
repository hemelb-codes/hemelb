HEMELB_CODE = ../../../../../Code

CXXFLAGS = -std=c++11 -g -I$(HEMELB_CODE)
LDLIBS = -lcppunit

objects = Vector.o Vector3D.o TriangleSorter.o test.o

test : $(objects)
	$(CXX) $(LDLIBS) $(objects) -o $@

test.o : test.cpp test_*.hpp

Vector.o : Vector.cpp

Vector3D.o : $(HEMELB_CODE)/util/Vector3D.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

TriangleSorter.o : TriangleSorter.cpp
