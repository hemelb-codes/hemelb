
include config.mk

HEMELB_CODE = ../../../../../Code

CXXFLAGS = -std=c++11 -g -I$(HEMELB_CODE) -I$(BOOST_INCLUDE_DIR) -I$(CGAL_INCLUDE_DIR) -I$(HDF5_INCLUDE_DIR) -I$(VTK_INCLUDE_DIR) 

vtk_needed_libs = vtkIOXML vtkCommonDataModel vtkCommonCore
vtk_libs = $(foreach vlib,$(vtk_needed_libs),-l$(vlib)-$(VTK_VERSION))

LDLIBS =  -lcppunit -L$(BOOST_LIB_DIR) -lboost_system -lboost_filesystem  -L$(CGAL_LIBRARY_DIR) -lCGAL -lgmp -L$(HDF5_LIB_DIR) -lhdf5 -L$(VTK_LIBRARY_DIR) -Wl,-rpath,$(VTK_LIBRARY_DIR) $(vtk_libs)  

lib_objects = Vector.o Vector3D.o TriTree.o TriangleSorter.o MkCgalMesh.o \
	Neighbours.o SegmentFactory.o SurfaceVoxeliser.o FloodFill.o \
	SectionTree.o SectionTreeBuilder.o H5.o

all : test voxer_bench

test : $(lib_objects) test.o
	$(CXX) $(LDLIBS) $^ -o $@

voxer_bench : $(lib_objects) voxer_bench.o
	$(CXX) $(LDLIBS) $^ -o $@

test.o : test.cpp test_*.hpp TestResources/*

Vector.o : Vector.cpp

Vector3D.o : $(HEMELB_CODE)/util/Vector3D.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

TriangleSorter.o : TriangleSorter.cpp

.PHONY : clean
clean :
	-rm test test.o voxer_bench voxer_bench.o $(lib_objects)
