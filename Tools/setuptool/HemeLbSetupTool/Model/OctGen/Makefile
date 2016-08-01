
include config.mk

HEMELB_CODE = ../../../../../Code

CXXFLAGS = -std=c++11 -g -I$(HEMELB_CODE) -I$(VTK_INCLUDE_DIR)
vtk_needed_libs = vtkIOXML vtkCommonDataModel vtkCommonCore
vtk_libs = $(foreach vlib,$(vtk_needed_libs),-l$(vlib)-$(VTK_VERSION))

LDLIBS = -L$(VTK_LIBRARY_DIR) -Wl,-rpath,$(VTK_LIBRARY_DIR) $(vtk_libs) -lcppunit

objects = Vector.o Vector3D.o TriTree.o TriangleSorter.o SurfaceVoxeliser.o \
	Neighbours.o VoxelClassifier.o test.o

test : $(objects)
	$(CXX) $(LDLIBS) $(objects) -o $@

test.o : test.cpp test_*.hpp

Vector.o : Vector.cpp

Vector3D.o : $(HEMELB_CODE)/util/Vector3D.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

TriangleSorter.o : TriangleSorter.cpp

.PHONY : clean
clean :
	-rm test $(objects)
